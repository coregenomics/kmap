library(magrittr)
library(Biostrings)
library(GenomicRanges)
library(BiocParallel)

## Allow BSgenome to be converted to GRanges object.
setAs("BSgenome", "GRanges",
      function(from) {
          makeGRangesFromDataFrame(
              data.frame(chr = seqnames(from),
                         start = rep(1, length(from)),
                         end = seqlengths(from)),
              ignore.strand = TRUE) %>%
              `names<-`(NULL)
      })

## Allow BSgenome to be converted to Views via GRanges.
setAs("BSgenome", "Views",
      function(from) as(from, "GRanges") %>% BSgenomeViews(from, .))

## Get standard DNA bases from BSgenomeViews by masking off the
## non-standard bases and returning the resulting ranges.
gr_masked <- function(views, motif = "N") {
    views %>% as("DNAStringSet") %>% sapply(mask, motif) %>%
        as("RangesList") %>% shift(start(views) - 1) %>%
        `names<-`(seqnames(views)) %>%
        GRanges(seqinfo = seqinfo(views))
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

## FIXME yes, this function will be cleaned up.
mappable <- function() {
    library(BSgenome)
    bsgenome <- getBSgenome("BSgenome.Ecoli.NCBI.20080805")
    ## ## Sanity check - there should be some non-DNA bases
    ## alphabetFrequency(as(bsgenome, "Views"), baseOnly = TRUE)
    system.time(views <- stddna(bsgenome))
    system.time(kmers <- kmerize(views))
    ## ## Sanity check - there should be no non-standard bases,
    ## ## otherwise pdict creation will fail.
    ## alphabetFrequency(kmers, baseOnly = TRUE) %>% colSums()
    ## ## There isn't really much one can do to speed up the PDict creation.
    system.time(pdict <- PDict(as(kmers, "DNAStringSet")))
    ## This step is memory intensive.  Using even more than 1 CPU eats
    ## all the RAM!  Might need to chop up the DNAStringSet into
    ## smaller pieces.
    system.time(hits <- lapply(as(views, "DNAStringSet"), matchPDict,
                               pdict = pdict))
}
