suppressPackageStartupMessages({
    library(BSgenome)
    library(magrittr)
})

context("Masking non-standard DNA bases")

## Fixtures
##
## Biostring genome.  Generating a mock BSgenome relies on having many
## disk files, which is not a very practical fixture.  Therefore use
## the smallest existing genome with unknown bases.  This Ecoli has 3
## chromosomes with non-standard bases: NC008563, NC_004431, NC_002655
bsgenome <- getBSgenome("BSgenome.Ecoli.NCBI.20080805")
## BSgenomeView.  Using a simpler DNAString fixture is limited
## because it does not have a seqnames() accessor.
.seqname <- "NC_002655"
.width <- 20
.region_std <- 1
.get_start <- function(pattern, width = 3) {
    mask(bsgenome[[.seqname]], pattern = pattern) %>% gaps() %>%
        .[width(.) == width] %>% head(1) %>% start()
}
.region_start <- .get_start("N")
.region_middle <- .get_start("D", 1) - .width / 2
.region_end <- .get_start("N", 2) - .width + 2
.gr <- GRanges(seqnames = .seqname,
               ranges = IRanges(
                   start = c(.region_std, .region_start,
                             .region_middle, .region_end),
                   width = .width),
               seqinfo = seqinfo(bsgenome))
views <- Views(bsgenome, .gr)
views
## GRanges corresponding to view with no non-standard bases.
.dnas <- DNAStringSet(views) %>% as.character()
## Slow, but reliable way of checking for standard bases.
.dna_ir <- function(dna_str) {
    bases <- unlist(strsplit(dna_str, split = NULL))
    bases %in% DNA_BASES %>% Rle() %>% IRanges()
}
gr <- lapply(.dnas, .dna_ir) %>%
    IRangesList %>%
    shift(start(.gr) - 1) %>%
    `names<-`(seqnames(.gr)) %>%
    GRanges(seqinfo = seqinfo(bsgenome))

## Function to test:
##   stddna_from_views(bsgenomeviews)
## Input:
##   BSgenomeViews
## Returns:
##   GRanges of regions containing only standard DNA bases (namely,
##   the DNA_BASES variable in the BSgenome packages)
## Description:
##   Create IRanges of standard DNA bases (A, C, G, T).
test_that("stddna returns GRanges-class for Views", {
    expect_is(stddna_from_views(views), "GRanges")
})
test_that("stddna returns empty XRanges for empty Views", {
    gr_ <- GRanges(seqinfo = seqinfo(bsgenome))
    views_ <- Views(bsgenome, gr_)
    expect_is(stddna_from_views(views_), "GRanges")
    expect_equal(stddna_from_views(views_), gr_)
})
test_that("stddna returns correct GRanges value for Views", {
    expect_equal(sort(stddna_from_views(views)),
                 sort(gr))
})
test_that("stddna returns contiguous GRanges", {
    views_ <- views[c(1, length(views))]
    gr_ <- gr[c(1, length(gr))]
    expect_equal(sort(stddna_from_views(views_)),
                 sort(gr_))
})
