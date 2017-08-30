context("Masking non-standard DNA bases")

test_that("stddna returns GRanges-class for Views", {
    expect_is(stddna_from_views(views), "GRanges")
})
test_that("stddna validates bsgenome", {
    expect_error(stddna_from_genome(NULL), "BSgenomeViews")
    expect_error(stddna_from_genome(GRanges()), "BSgenomeViews")
    expect_error(stddna_from_genome("sacCer2"), "BSgenomeViews")
    expect_error(stddna_from_genome("BSgenome.Scerevisiae.UCSC.sacCer2"),
                 "BSgenomeViews")
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
