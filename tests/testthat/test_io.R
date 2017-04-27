context("Generation of alignment files")
source("fixtures.R")

## Function to test:
##   save_as_sample_file(filename, views)
test_that("save_as_sample_file is nearly serializable", {
    ## We can't serialize genomic information, but can validate.
    filename <- tempfile()
    on.exit(unlink(filename))
    save_as_sample_file(filename, views)
    df <- read.delim(filename, header = FALSE)
    gr <- GRanges(df[, 1])
    dna <- df[, 2] %>% as.character()
    expect_equivalent(gr, granges(views))
    expect_equivalent(dna, views %>% as.character())
})
