context("Utils")

test_that("ir2gr ensures inputs have the same length", {
    ## in a section of the Ecoli genome.
    ir <- ranges(views)
    expect_error(ir2gr(head(ir, 2), views))
})
