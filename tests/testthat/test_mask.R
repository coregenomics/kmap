context("Masking non-standard DNA bases")

## Function to test:
##   stddna_from_views(bsgenomeviews)
## Input:
##   BSgenomeViews
## Returns:
##   GRanges of regions containing only standard DNA bases (namely,
##   the DNA_BASES variable in the BSgenome packages)
## Description:
##   Create GRanges of standard DNA bases (A, C, G, T).
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
