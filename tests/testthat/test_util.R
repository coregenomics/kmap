context("Utils")

test_that("ir2gr ensures inputs have the same length", {
    ## in a section of the Ecoli genome.
    ir <- ranges(views)
    expect_error(ir2gr(head(ir, 2), views))
})

test_that("coerce from BSgenome returns GRanges", {
    result <- as(bsgenome, "GRanges")
    expect_s4_class(result, "GRanges")
})

test_that("coerce from BSgenome returns BSgenomeViews", {
    result <- as(bsgenome, "BSgenomeViews")
    expect_s4_class(result, "BSgenomeViews")
})
