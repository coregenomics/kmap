context("Test mappable")

## Fixtures
genome <- "BSgenome.Scerevisiae.UCSC.sacCer2"
path_cache <- tempdir()
clear_cache <- function() {
    unlink(tempdir(), recursive = TRUE)
    dir.create(tempdir())
}

test_that("mappable validates genome name", {
    expect_error(mappable(genome = "nonsense"), "BSgenome")
    expect_error(mappable(genome = NULL), "BSgenome")
})

test_that("mappable loads GRanges from the the file cache", {
    ## Start with an empty cache in case of any failed or interrupted tests.
    clear_cache()
    on.exit(clear_cache())
    ## Write test data to the cache.
    genome_short <- genome %>% strsplit("\\.") %>% unlist() %>% tail(1)
    gr <- GRanges("chrI:1000-1200", seqinfo = Seqinfo(genome = genome_short))
    name <- paste("kmap", formals(mappable)$kmer, genome, sep = "_")
    bfc <- BiocFileCache(path_cache)
    suppressWarnings(path_gff3_1 <- bfcnew(bfc, name, ext = "gff3"))
    rtracklayer::export(gr, path_gff3_1)
    ## Read from cache.
    gr_from_cache <- mappable(genome, cache_path = path_cache)
    mcols(gr_from_cache) <- NULL
    expect_equal(gr_from_cache, gr)

    ## Multiple files of the same name in the cache should raise an error.
    suppressWarnings(path_gff3_2 <- bfcnew(bfc, name, ext = "gff3"))
    rtracklayer::export(gr, path_gff3_2)
    expect_error(mappable(genome, cache_path = path_cache), "more than 1")
})

test_that("mappable generates and returns consistent GRanges", {
    ## Travis seems to kill this test; perhaps for taking too long?
    skip_on_travis()
    ## TODO: Add ability to use fasta file for genome to speed up this test.
    clear_cache()
    on.exit(clear_cache())
    expect_message(gr_1 <- mappable(genome, cache_path = path_cache),
                   "Saving result")
    expect_s4_class(gr_1, "GRanges")
    expect_message(gr_2 <- mappable(genome, cache_path = path_cache),
                   "Reading the cached")
    expect_s4_class(gr_2, "GRanges")
    mcols(gr_2) <- NULL
    expect_equal(gr_2, gr_1)
})
