## Biostring genome.  Generating a mock BSgenome relies on having many
## disk files, which is not a very practical fixture.  Therefore use
## the smallest existing genome with unknown bases.  This Ecoli has 3
## chromosomes with non-standard bases: NC008563, NC_004431, NC_002655
bsgenome <- BSgenome::getBSgenome("BSgenome.Ecoli.NCBI.20080805")

## BSgenomeView.  Using a simpler DNAString fixture is limited
## because it does not have a seqnames() accessor.
.seqname <- "NC_002655"
.width <- 20
.region_std <- 1
.get_start <- function(pattern, width = 3) {
    mask(bsgenome[[.seqname]], pattern = pattern) %>% gaps() %>%
        .[width(.) == width] %>% head(1) %>% start()
}
.region_start <- .get_start("N")
.region_middle <- .get_start("D", 1) - .width / 2
.region_end <- .get_start("N", 2) - .width + 2
.gr <- GRanges(seqnames = .seqname,
               ranges = IRanges(
                   start = c(.region_std, .region_start,
                             .region_middle, .region_end),
                   width = .width),
               seqinfo = seqinfo(bsgenome))
views <- Views(bsgenome, .gr)
views

## GRanges corresponding to view with no non-standard bases.
.dnas <- DNAStringSet(views) %>% as.character()
## Slow, but reliable way of checking for standard bases.
.dna_ir <- function(dna_str) {
    bases <- unlist(strsplit(dna_str, split = NULL))
    bases %in% DNA_BASES %>% Rle() %>% IRanges()
}
gr <- lapply(.dnas, .dna_ir) %>%
    IRangesList %>%
    shift(start(.gr) - 1) %>%
    `names<-`(seqnames(.gr)) %>%
    GRanges(seqinfo = seqinfo(bsgenome))

## kmer size smaller than gr .width
kmer <- 10
