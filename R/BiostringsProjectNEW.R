# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("GenomicRanges")
# biocLite("Biostrings")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# # biocLite("BiocParallel")
# biocLite("BSgenome.Celegans.UCSC.ce2")


library(Biostrings)
library(GenomicRanges)
# library(BiocParallel)
# library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
# library(BiocInstaller)
# biocLite("BSgenome.Ecoli.NCBI.20080805")
# biocLite("BSgenome.Scerevisiae.UCSC.sacCer2")

# # Subset to chromosome x for quick computation


# Create blacklist of unmappable regions, mappable regions w/ only ACGT

# # Determine starts for reads of width 36, create views

stddna_chrom <- function(dna_str){dna_ir <- dna_str %>% mask("N") %>% as("IRanges")}
stddna <- function(bsgenome){
    ranges <- makeGRangesFromDataFrame(
        data.frame(chr = seqnames(bsgenome),
                   start = rep(1, length(bsgenome)),
                   end = seqlengths(bsgenome)),ignore.strand = TRUE) %>%
        `names<-`(NULL)
    seqs <- Views(bsgenome, ranges) %>% as("DNAStringSet")
    genome_ir <- sapply(seqs, stddna_chrom) %>% RangesList %>%
        `names<-`(seqnames(bsgenome))%>%GRanges
}
