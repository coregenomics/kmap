context("Generating query dictionary for mapping")

## Function to test:
##   kmerize(views, kmer)
test_that("kmerize returns BSgenomeViews-class", {
     expect_is(kmerize(views, kmer), "BSgenomeViews")
})
test_that("kmerize returns ranges no wider than the kmer size", {
    kmer_ <- 10
    expect_lte(width(kmerize(views, kmer_)) %>% max(), kmer_)
    kmer_ <- 100
    expect_lte(width(kmerize(views, kmer_)) %>% length() %>% max(), kmer_)
    kmer_ <- 1
    expect_lte(width(kmerize(views, kmer_)) %>% max(), kmer_)
})
test_that("kmerize throws error for invalid kmer size", {
    expect_error(kmerize(views, -1), "kmer")
    expect_error(kmerize(views, 0), "kmer")
    expect_error(kmerize(views, 10.5), "kmer")
})
