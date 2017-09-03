## Biostring genome.  Generating a mock BSgenome relies on having many
## disk files, which is not a very practical fixture.  Therefore use
## the smallest existing genome with unknown bases.  This Ecoli has 3
## chromosomes with non standard bases: NC008563, NC_004431, NC_002655
genome <- "BSgenome.Ecoli.NCBI.20080805"
bsgenome <- BSgenome::getBSgenome(genome)
## BSgenomeView.  Using a simpler DNAString fixture is limited
## because it does not have a seqnames accessor.
.seqname <- "NC_002655"
.width <- 20
.region_std <- 1
.get_start <- function(pattern, width = 3, bsgenome_ = bsgenome,
                       .seqname_ = .seqname) {
    xstringviews <- Biostrings::mask(bsgenome_[[.seqname_]], pattern = pattern)
    xstringviews <- GenomicRanges::gaps(xstringviews)
    xstringviews[GenomicRanges::width(xstringviews) == width]
    GenomicRanges::start(head(xstringviews, 1))
}
.region_start <- .get_start("N")
.region_middle <- .get_start("D", 1) - .width / 2
.region_end <- .get_start("N", 2) - .width + 2
.gr <- GenomicRanges::GRanges(seqnames = .seqname,
                              ranges = IRanges::IRanges(
                                  start = c(.region_std, .region_start,
                                            .region_middle, .region_end),
                                  width = .width),
                              seqinfo = GenomeInfoDb::seqinfo(bsgenome))
views <- BSgenome::Views(bsgenome, .gr)
views

## GRanges corresponding to view with no non standard bases.
.dnas <- as.character(Biostrings::DNAStringSet(views))
## Slow, but reliable way of checking for standard bases.
.dna_ir <- function(dna_str) {
    bases <- unlist(strsplit(dna_str, split = NULL))
    bases_std <- bases %in% Biostrings::DNA_BASES
    bases_std <- S4Vectors::Rle(bases_std)
    IRanges::IRanges(bases_std)
}
gr <- IRanges::IRangesList(lapply(.dnas, .dna_ir))
gr <- GenomicRanges::shift(gr, start(.gr) - 1)
names(gr) <- GenomeInfoDb::seqnames(.gr)
gr <- GenomicRanges::GRanges(gr, seqinfo = GenomeInfoDb::seqinfo(bsgenome))

## kmer size smaller than gr .width
kmer <- 10
