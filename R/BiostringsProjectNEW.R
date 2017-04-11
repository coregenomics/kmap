#' @importFrom BSgenome BSgenome BSgenomeViews subject elementNROWS
#' @importFrom BiocParallel bplapply
#' @importFrom Biostrings DNAStringSet DNA_ALPHABET DNA_BASES PDict matchPDict mask
#' @importFrom GenomeInfoDb seqinfo seqlengths
#' @importFrom GenomicRanges GRanges end gaps reduce seqnames shift start width
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors List runLength<-
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @importFrom plyr .

## Allow BSgenome to be converted to GRanges object.
setAs("BSgenome", "GRanges",
      function(from) {
          GenomicRanges::makeGRangesFromDataFrame(
              data.frame(chr = seqnames(from),
                         start = rep(1, length(from)),
                         end = seqlengths(from)),
              ignore.strand = TRUE) %>%
              `names<-`(NULL)
      })

## Allow BSgenome to be converted to Views via GRanges.
setAs("BSgenome", "Views",
      function(from) as(from, "GRanges") %>% BSgenomeViews(from, .))

## Convert IRanges back to GRanges using metadata from Views.
ir2gr <- function(ir, views) {
    ir %>% shift(start(views) - 1) %>%
         `names<-`(seqnames(views)) %>%
         GRanges(seqinfo = seqinfo(views))
}

## Get standard DNA bases from BSgenomeViews by masking off the
## non-standard bases and returning the resulting ranges.
gr_masked <- function(views, motif = "N") {
    views %>% as("DNAStringSet") %>% sapply(mask, motif) %>%
        as("RangesList") %>% ir2gr(views)
}

stddna_from_views <- function(views) {
    bases_nonstd <- setdiff(DNA_ALPHABET, DNA_BASES)
    ## FIXME maybe parallelize this for larger genomes.
    bplapply(bases_nonstd, gr_masked, views = views) %>% sapply(gaps) %>%
        List() %>% unlist() %>% reduce() %>% gaps()
}

## BSgenomeViews of standard DNA from BSgenome object.
stddna_from_genome <- function(bsgenome) {
    ## Convert to Views here for convenience, instead of subsetting by
    ## chromosomes.  Views allows more vector operations and preserves
    ## metadata better than shuttling around chromosomes.
    bsgenome %>% as("Views") %>% stddna_from_views() %>%
        BSgenomeViews(subject = bsgenome)
}

## Chop up BSgenomeViews into kmers.
kmerize <- function(views, kmer = 36) {
    ## Ensure that end(views) >= start(views).  Drop rows where this
    ## is not the case.
    views <- views[width(views) > kmer]
    ## Optimize runtime by disabling USE.NAMES.
    starts <- mapply(seq,
                     start(views),
                     end(views) - kmer,
                     USE.NAMES = FALSE)
    gr_lengths <- seqnames(views) %>%
        as.data.frame() %>%
        dplyr::mutate(len = sapply(starts, length)) %>%
        dplyr::group_by(value) %>%
        dplyr::summarise(lengths = sum(len)) %>%
        .$lengths
    gr_seqnames <- seqnames(views) %>% `runLength<-`(gr_lengths)
    gr <- GRanges(ranges = IRanges(start = starts %>% unlist(),
                                   width = kmer),
                  seqnames = gr_seqnames)
    BSgenomeViews(subject(views), gr)
}

.range_hits <- function(xstring, pdict) {
    matches <- matchPDict(pdict, xstring)
    hits <- elementNROWS(matches) > 1
    matches %>% as("CompressedIRangesList") %>% .[hits] %>%
        unlist() %>% IRanges::reduce()
}

ranges_hits <- function(views, pdict, indices = NULL) {
    if (is.null(indices)) indices <- 1:length(views)
    views_sub <- views[indices]
    irl <- lapply(views_sub, .range_hits, pdict) %>% List()
    ## matchPDict() converted the BSgenomeViews into DNAStringSets,
    ## and therefore the start offsets and seqinfo was lost.  We need
    ## to readd that information using ir2gr() before combining all
    ## the ranges.
    lapply(seq_along(irl), function(i) ir2gr(irl[i], views_sub[i])) %>%
        List() %>% unlist() %>% GenomicRanges::reduce()
}

## FIXME yes, this function will be cleaned up.
mappable <- function() {
    message(sprintf(
        "[%6.2f sec] to load the genome",
        system.time(
            bsgenome <-
                BSgenome::getBSgenome("BSgenome.Ecoli.NCBI.20080805"))[3]))
    ## ## Sanity check - there should be some non-DNA bases
    ## alphabetFrequency(as(bsgenome, "Views"), baseOnly = TRUE)
    message(sprintf(
        "[%6.2f sec] to subset to standard DNA bases in the genome",
        system.time(views <- stddna_from_genome(bsgenome))[3]))
    message(sprintf(
        "[%6.2f sec] to chop up the genome into 36-mers",
        system.time(kmers <- kmerize(views))[3]))
    ## ## Sanity check - there should be no non-standard bases,
    ## ## otherwise pdict creation will fail.
    ## alphabetFrequency(kmers, baseOnly = TRUE) %>% colSums()
    ## ## There isn't really much one can do to speed up the PDict creation.
    message(sprintf(
        "[%6.2f sec] to create the search dictionary",
        system.time(pdict <- PDict(as(kmers,
                                      "DNAStringSet")))[3]))
    ## Chop up the views into smaller pieces.
    message(sprintf(
        "[%6.2f sec] to search the genome",
        system.time(gr <- bplapply(seq_along(views),
                                   ranges_hits,
                                   views = views,
                                   pdict = pdict) %>%
                        List() %>%
                        unlist() %>%
                        GenomicRanges::reduce())[3]))
    gr
}
