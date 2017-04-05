context("Generating query dictionary for mapping")
source("fixtures.R")

## Function to test:
##   kmerize(views, kmer)
## Input:
##   BSgenomeViews
##   Integer
## Returns:
##   Biostrings PDict object to run sequence queries.
## Description:
##   Chop up genome into fixed size to run sequence queries.
## test_that("kmerize returns Views-class", {
##     expect_is(kmerize(views, kmer), "Views")
## })
test_that("kmerize returns ranges no wider than the kmer size", {
    skip("Not completed yet")
    result <- function() width(kmerize(views, kmer_))
    result_vec <- function() rep(kmer_ + 1, length(result_))
    kmer_ <- kmer
    expect_lt(result(), result_vec())
    kmer_ <- 100
    expect_lt(result(), result_vec())
    kmer_ <- 1
    expect_lt(result(), result_vec())
})

## Test negative value gives error
## Test zero value gives error
## Test floating point gives error
## Test 1-2 kmer sizes
## Test option to trim to constant size
