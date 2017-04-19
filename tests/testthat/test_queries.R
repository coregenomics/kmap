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
    kmer_ <- 10
    expect_lte(width(kmerize(views, kmer_)) %>% max(), kmer_)
    kmer_ <- 100
    expect_lte(width(kmerize(views, kmer_)) %>% length() %>% max(), kmer_)
    kmer_ <- 1
    expect_lte(width(kmerize(views, kmer_)) %>% max(), kmer_)
})

## Test negative value gives error
## Test zero value gives error
## Test floating point gives error
## Test 1-2 kmer sizes
## Test option to trim to constant size
