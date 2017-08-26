context("Alignment")

## Function to test:
##   align(views, genome)
test_that("align generates .fasta files", {
    ## We can't use the NCBI BSgenome name directly as the QuasR::qAlign genome
    ## parameter, because the NCBI metadata has NAs in the publication date
    ## which confuses QuasR.  As a workaround, resave the genome as a file.
    file_genome <- tempfile(pattern = "genome-", fileext = ".fasta")
    on.exit(unlink(file_genome))
    ## Subset the genome using `head()` to speed up the calculation from several
    ## minutes to tens of seconds.
    writeXStringSet(as(bsgenome, "Views") %>% as("XStringSet") %>%
                    setNames(names(bsgenome)) %>% endoapply(head, 100),
                    file_genome)
    ## None of the `views` will map uniquely, so add in a uniquely mapping
    ## string to validate.
    gr_unique <- GRanges("NC_008563:1-100", seqinfo = seqinfo(views))
    views_1_unique <- Views(bsgenome, c(.gr, gr_unique))

    ## Rbowtie runs a subprocess, so one cannot suppress output in conventional
    ## ways.
    result <- align(views_1_unique, file_genome)
    expect_equal(result, gr_unique)
})
