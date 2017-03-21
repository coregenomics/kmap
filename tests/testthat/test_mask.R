suppressPackageStartupMessages({
    library(BSgenome)
    library(magrittr)
})

context("Masking non-standard DNA bases")

## Fixtures
##
## DNA string
.str <- "TTGAANNNAACTCNACTG"
dna_str <- DNAString(.str)
.bases <- unlist(strsplit(.str, split = NULL))
dna_ir <- IRanges(Rle(.bases %in% DNA_BASES))

## Generating a mock BSgenome relies on having many disk files, which
## is not a very practical fixture.  Therefore use the smallest
## existing genome with unknown bases.
ir_nonstd <- function(dnastring) {
    dnastring %>% maskMotif("N") %>% masks() %>% .[[1]]
}
count_nonstd <- function(dnastring) {
    ir_nonstd(dnastring) %>% width() %>% sum()
}
has_nonstd <- function(dnastring) {
    count_nonstd(dnastring) > 0
}
inspect_genome <- function(bsgenome) {
    for (chrom in seqnames(bsgenome)) {
        count <-  count_nonstd(bsgenome[[chrom]])
        print(c(chrom, count))
    }
}
## Ecoli has 3 chromosomes with non-standard bases: NC008563,
## NC_004431, NC_002655
bsgenome <- getBSgenome("BSgenome.Ecoli.NCBI.20080805")
##inspect_genome(bsgenome)
gr <- makeGRangesFromDataFrame(
    data.frame(chr = seqnames(bsgenome),
               start = rep(1, length(bsgenome)),
               end = seqlengths(bsgenome)),
    ignore.strand = TRUE) %>%
    `names<-`(NULL)

## Function to test:
##   stddna_chrom(x)
## Input:
##   DNAString
## Returns:
##   IRanges of regions containing only standard DNA bases (namely,
##   the DNA_BASES variable in the BSgenome packages)
## Description:
##   Create IRanges of standard DNA bases (A, C, G, T).
test_that("stddna returns IRanges-class for DNAString", {
    expect_is(stddna_chrom(dna_str), "IRanges")
})
test_that("stddna returns empty IRanges for empty DNAString", {
    dna_str <- DNAString()
    dna_ir <- IRanges()
    expect_is(stddna_chrom(dna_str), "IRanges")
    expect_equal(stddna_chrom(dna_str), dna_ir)
})
test_that("stddna returns correct IRanges value for DNAString", {
    expect_equal(stddna_chrom(dna_str), dna_ir)
})

## Function to test:
##   stddna(genome)
## Input:
##   BSgenome
## Returns:
##   GRanges of regions containing only standard DNA bases (namely,
##   the DNA_BASES variable in the BSgenome packages)
## Description:
##   Create GRanges of standard DNA bases (A, C, G, T).  At the time
##   of writing, BSgenome does not allow subsetting DNAString objects
##   of more than one chromosome at a time, therefore this function
##   calls stddna() for each chromosome.  This should be a fast
##   operation and should not need parallel library computation.
test_that("stddna returns GRanges-class for BSgenome", {
    expect_is(stddna(bsgenome), "GRanges")
})
## It's not practical to create an empty BSgenome, as the BSgenome
## object requires disk files.  Therefore test a completely sequenced
## genome with no non-standard DNA bases.
test_that("stddna returns contiguous GRanges for fully sequenced BSgenome", {
    bsgenome <- getBSgenome("BSgenome.Scerevisiae.UCSC.sacCer2")
    gr <- makeGRangesFromDataFrame(
        data.frame(chr = seqnames(bsgenome),
                   start = rep(1, length(bsgenome)),
                   end = seqlengths(bsgenome)),
        ignore.strand = TRUE) %>%
        `names<-`(NULL)
    expect_is(stddna(bsgenome), "GRanges")
    expect_equal(stddna(bsgenome), gr)
})
test_that("stddna returns correct IRanges value for BSgenome", {
    expect_equal(stddna(bsgenome), gr)
})
