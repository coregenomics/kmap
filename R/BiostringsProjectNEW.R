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
# genes <- c(BSgenome.Hsapiens.UCSC.hg19[[19]],BSgenome.Hsapiens.UCSC.hg19[[20]])

# Create blacklist of unmappable regions, mappable regions w/ only ACGT
# unmappable <- genes %>% mask("N") %>% gaps %>% as("IRanges")
# mappable <- genes %>% mask("N") %>% as("IRanges")
#
# # Determine starts for reads of width 36, create views
# starts_map <- mapply(seq,start(mappable),end(mappable)-36)
# views_map <- IRanges::Views(genes, start=unlist(starts_map), width=36)
# queries <- DNAStringSet(views_map)
# system.time({PDict1 <- PDict(head(queries,100000))
# alignment <- matchPDict(PDict1, chr19, max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, fixed = TRUE, algorithm = "auto", verbose = FALSE)
# hits <- sapply(alignment,length)>1
# hit_sequences <- as(alignment, "CompressedIRangesList")%>%.[hits]%>%unlist%>%IRanges::reduce()
# })
# hit_sequences
#
# matches <- function(query){PDict1 <- PDict(query)
# alignment <- matchPDict(PDict1, genes, max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, fixed = TRUE, algorithm = "auto", verbose = FALSE)
# hits <- sapply(alignment,length)>1
# hit_sequences <- as(alignment, "CompressedIRangesList")%>%.[hits]%>%unlist%>%IRanges::reduce()
# }
# queries_list <- split(queries, ceiling(seq_along(queries)/100000))
# MulticoreParam(progressbar = TRUE)
# hits <- bplapply(head(queries_list), matches)
# unmappable_ranges <- hits%>%c(unmappable)%>%as("List")%>%unlist()%>%IRanges::reduce()

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
