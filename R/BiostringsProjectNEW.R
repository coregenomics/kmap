#' @importFrom BSgenome BSgenome BSgenomeViews subject elementNROWS
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom Biostrings DNAStringSet DNA_ALPHABET DNA_BASES PDict matchPDict mask
#' @importFrom GenomeInfoDb seqinfo seqlengths
#' @importFrom GenomicRanges GRanges end gaps reduce seqnames shift start width
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors List runLength<-
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @importFrom plyr .

## The `setAs()` functions are intentionally not documented by themselves so as
## to not clobber the base R `as()` help.  Bioconductor core packages typically
## document `as()` alongside the class descriptions / constructors.  Here we're
## adding `BSgenome` conversions which might not be of interest to upstream.

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

#' Subset to standard DNA bases
#'
#' Functions for subsetting to standard DNA for a BSgenome (or BSgenomeViews).
#'
#' @param bsgenome The \code{\link[BSgenome]{BSgenome}} DNA to subset.
#' @param views The \code{\link[BSgenome]{BSgenomeViews}} DNA to subset.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation,
#'     or a \code{\link[base]{list}} of
#'     \code{\link[BiocParallel]{BiocParallelParam}} instances, to be applied in
#'     sequence for nested calls to \code{BiocParallel} functions.
#' @return A \code{BSgenomeViews} object with non-standard DNA bases (i.e. bases
#'     not in \code{\link[BSgenome]{DNA_BASES}}) removed for
#'     \code{\link{stddna_from_genome}}, and a \code{GRanges} object for
#'     \code{\link{stddna_from_views}}.
#' @name stddna
NULL

#' @rdname stddna
stddna_from_genome <- function(bsgenome, BPPARAM = bpparam()) {
    ## Convert to Views here for convenience, instead of subsetting by
    ## chromosomes.  Views allow more vector operations and preserve metadata
    ## better than shuttling around chromosomes.
    bsgenome %>% as("Views") %>% stddna_from_views(BPPARAM = BPPARAM) %>%
        BSgenomeViews(subject = bsgenome)
}

#' @rdname stddna
stddna_from_views <- function(views, BPPARAM = bpparam()) {
    bases_nonstd <- setdiff(DNA_ALPHABET, DNA_BASES)
    ## FIXME maybe parallelize this for larger genomes.
    bplapply(bases_nonstd, gr_masked, views = views, BPPARAM = BPPARAM) %>%
        sapply(gaps) %>% List() %>% unlist() %>% reduce() %>% gaps()
}

#' Remove motifs from BSgenomeViews
#'
#' Filter out DNA motifs to create a \code{\link[GenomicRanges]{GRanges}}
#' instance.
#'
#' @inheritParams stddna
#' @param motif The motif to remove in the sequences.
#' @return A \code{\link[GenomicRanges]{GRanges}} object with \code{motif}
#'     sequences removed.
gr_masked <- function(views, motif = "N") {
    views %>% as("DNAStringSet") %>% sapply(mask, motif) %>%
        as("RangesList") %>% ir2gr(views)
}

#' Convert IRanges back to GRanges using metadata from Views.
#'
#' We lose the \code{\link[GenomicRanges]{seqnames}} and
#' \code{\link[GenomeInfoDb]{seqinfo}} metadata when transforming the underlying
#' \code{\link[Biostrings]{DNAString}} of a
#' \code{\link[BSgenome]{BSgenomeViews}} object into
#' \code{\link[IRanges]{IRanges}}.  \code{\link{ir2gr}} recovers that
#' information from the original \code{\link[BSgenome]{BSgenomeViews}} to up
#' convert \code{\link[IRanges]{IRanges}} to
#' \code{\link[GenomicRanges]{GRanges}}.
#'
#' @inheritParams stddna
#' @param ranges IRanges
#' @return GRanges
#' @examples
#' bsgenome <- BSgenome::getBSgenome("BSgenome.Ecoli.NCBI.20080805")
#' views <- BSgenomeViews(bsgenome, GRanges("AC_000091:11-20"))
#' views
#' remaining <- mask(views[[1]], "T")
#' ir <- as(remaining, "IRanges")
#' ir
#' ## Convert IRanges to GRanges
#' ir2gr(ir, views)
ir2gr <- function(ranges, views) {
    ranges %>% shift(start(views) - 1) %>%
         `names<-`(seqnames(views)) %>%
         GRanges(seqinfo = seqinfo(views))
}

#' Chop up BSgenomeViews into overlapping k-mers.
#'
#' @param views The \code{\link[BSgenome]{BSgenomeViews}} sequences to segment
#'     into overlapping k-mers.
#' @return A \code{\link[BSgenome]{BSgenomeViews}} object with
#'     \code{\link[GenomicRanges]{width}} exactly equal to \code{kmer} size.
#'     Ranges smaller than the \code{kmer} are dropped.
#' @export
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

range_hits <- function(xstring, pdict) {
    matches <- matchPDict(pdict, xstring)
    hits <- elementNROWS(matches) > 1
    matches %>% as("CompressedIRangesList") %>% .[hits] %>%
        unlist() %>% IRanges::reduce()
}

ranges_hits <- function(views, pdict, indices = NULL) {
    if (is.null(indices)) indices <- 1:length(views)
    views_sub <- views[indices]
    irl <- lapply(views_sub, range_hits, pdict) %>% List()
    ## matchPDict() converted the BSgenomeViews into DNAStringSets, and
    ## therefore the start offsets and seqinfo was lost.  We need to readd that
    ## information using ir2gr() before combining all the ranges.
    lapply(seq_along(irl), function(i) ir2gr(irl[i], views_sub[i])) %>%
        List() %>% unlist() %>% GenomicRanges::reduce()
}

## FIXME yes, this function will be cleaned up.
#' @export
mappable <- function(genome, BPPARAM = bpparam()) {
    message(sprintf(
        "[%6.2f sec] to load the genome",
        system.time(
            bsgenome <-
                BSgenome::getBSgenome(genome))[3]))
    ## ## Sanity check - there should be some non-DNA bases
    ## alphabetFrequency(as(bsgenome, "Views"), baseOnly = TRUE)
    message(sprintf(
        "[%6.2f sec] to subset to standard DNA bases in the genome",
        system.time(views <- stddna_from_genome(bsgenome, BPPARAM))[3]))
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
                                   pdict = pdict,
                                   BPPARAM = BPPARAM) %>%
                        List() %>%
                        unlist() %>%
                        GenomicRanges::reduce())[3]))
    gr
}
