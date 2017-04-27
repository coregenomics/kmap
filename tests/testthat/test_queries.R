context("Generating query dictionary for mapping")
source("fixtures.R")

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

## Test negative value gives error
test_that("negative value of kmer size", {
    expect_error(kmerize(views, -1), "cannot have a negative kmer size")
})

## Test zero value gives error
test_that("zero value of kmer size", {
    expect_error(kmerize(views, 0), "cannot have a zero value of kmer")
})

## Test floating point gives error
## Test 1-2 kmer sizes
test_that("")
kmer_1 <- 19
kmer
## Test option to trim to constant size

