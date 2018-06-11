context("Utils")

test_that("expand_rle grows input Rle", {
    input <- Rle(c("a", "a", "a", "b", "c"))
    size <- Rle(c(1, 1, 2, 1, 3))
    expected <- Rle(c("a", "a", "a", "a", "b", "c", "c", "c"))
    result <- expand_rle(input, size)
    expect_equal(expected, result)
    ## Test case of duplicated run values b and c both being run length 2.
    input <- Rle(c("a", "a", "a", "b", "c", "c"))
    size <- Rle(c(1, 2, 1, 2, 2, 1))
    expected <- Rle(c("a", "a", "a", "a", "b", "b", "c", "c", "c"))
    result <- expand_rle(input, size)
    expect_equal(expected, result)
})

test_that("coerce from BSgenome returns GRanges", {
    result <- as(bsgenome, "GRanges")
    expect_s4_class(result, "GRanges")
})

test_that("coerce from BSgenome returns BSgenomeViews", {
    result <- as(bsgenome, "BSgenomeViews")
    expect_s4_class(result, "BSgenomeViews")
})
