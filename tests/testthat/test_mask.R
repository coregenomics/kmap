suppressPackageStartupMessages({
    library(BSgenome)
    library(magrittr)
})

context("Masking non-standard DNA bases")

## Functions to create Fixtures
##
## Allow BSgenome to be converted to GRanges object.
setAs("BSgenome", "GRanges",
      function(from) {
          makeGRangesFromDataFrame(
              data.frame(chr = seqnames(from),
                         start = rep(1, length(from)),
                         end = seqlengths(from)),
              ignore.strand = TRUE)
      })
## Allow BSgenome to be converted to DNAStringSet.
setAs("BSgenome", "Views",
      function(from) as(from, "GRanges") %>% BSgenomeViews(from, .))
## Find non-standard DNA bases.
ir_nonstd <- function(dnastring) {
    dnastring %>% maskMotif("N") %>% masks() %>% .[[1]] # nolint
}
## Find non-standard DNA bases.
gr_nonstd <- function(bsgenome) {
    ## Use BSgenomeViews To iterate over the chromosomes.
    lapply(as(bsgenome, "Views"), ir_nonstd) %>%
        as("RangesList") %>%
        as("GRanges")
}

## Fixtures
##
## DNA string
.str <- "TTGAANNNAACTCNACTG"
dna_str <- DNAString(.str)
.bases <- unlist(strsplit(.str, split = NULL))
dna_ir <- IRanges(Rle(.bases %in% DNA_BASES))
## Biostring genome
##
## Generating a mock BSgenome relies on having many disk files, which
## is not a very practical fixture.  Therefore use the smallest
## existing genome with unknown bases.  This Ecoli has 3 chromosomes
## with non-standard bases: NC008563, NC_004431, NC_002655
bsgenome <- getBSgenome("BSgenome.Ecoli.NCBI.20080805")
gr <- GenomicRanges::setdiff(as(bsgenome, "GRanges"),
                             gr_nonstd(bsgenome))

## Function to test:
##   stddna_chrom(bsgenomeviews)
## Input:
##   BSgenomeViews
## Returns:
##   IRanges of regions containing only standard DNA bases (namely,
##   the DNA_BASES variable in the BSgenome packages)
## Description:
##   Create IRanges of standard DNA bases (A, C, G, T).
test_that("stddna returns IRanges-class for Views", {
    expect_is(stddna_chrom(dna_str), "IRanges")
})
test_that("stddna returns empty IRanges for empty Views", {
    dna_str <- DNAString()
    dna_ir <- IRanges()
    expect_is(stddna_chrom(dna_str), "IRanges")
    expect_equal(stddna_chrom(dna_str), dna_ir)
})
test_that("stddna returns correct IRanges value for Views", {
    expect_equal(stddna_chrom(dna_str), dna_ir)
})

## Function to test:
##   stddna(genome)
## Input:
##   BSgenome
## Returns:
##   GRanges of regions containing only standard DNA bases (namely,
##   the DNA_BASES variable in the BSgenome packages)
## Description:
##   Create GRanges of standard DNA bases (A, C, G, T).  At the time
##   of writing, BSgenome does not allow subsetting DNAString objects
##   of more than one chromosome at a time, therefore this function
##   calls stddna() for each chromosome.  This should be a fast
##   operation and should not need parallel library computation.
test_that("stddna returns GRanges-class for BSgenome", {
    expect_is(stddna(bsgenome), "GRanges")
})
## It's not practical to create an empty BSgenome, as the BSgenome
## object requires disk files.  Therefore test a completely sequenced
## genome with no non-standard DNA bases.
test_that("stddna returns contiguous GRanges for fully sequenced BSgenome", {
    bsgenome <- getBSgenome("BSgenome.Scerevisiae.UCSC.sacCer2")
    gr <-GenomicRanges::setdiff(as(bsgenome, "GRanges"),
                                gr_nonstd(bsgenome))
    ## Make sure our fixture contains no non-standard DNA bases.
    expect_equal(gr, as(bsgenome, "GRanges") %>% `names<-`(NULL))
    ## Actual tests.
    expect_is(stddna(bsgenome), "GRanges")
    expect_equal(stddna(bsgenome), gr)
})
test_that("stddna returns correct IRanges value for BSgenome", {
    expect_equal(stddna(bsgenome), gr)
})
