#' @importFrom BSgenome BSgenome BSgenomeViews subject elementNROWS
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom Biostrings DNAStringSet DNA_ALPHABET DNA_BASES PDict
#'     writeXStringSet matchPDict mask
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb seqinfo seqlengths
#' @importFrom GenomicRanges GRanges end gaps reduce seqnames shift start width
#' @importFrom IRanges IRanges
#' @importFrom QuasR alignments qAlign
#' @importFrom S4Vectors List runLength<-
#' @importFrom magrittr %>%
#' @importFrom methods as
#' @importFrom plyr .

## The `setAs()` functions are intentionally not documented by themselves so as
## to not clobber the base R `as()` help.  Bioconductor core packages typically
## document `as()` alongside the class descriptions / constructors.  Here we're
## adding `BSgenome` conversions which might not be of interest to upstream.

## Allow BSgenome to be converted to GRanges object.
setAs("BSgenome", "GRanges",
      function(from) {
          GenomicRanges::makeGRangesFromDataFrame(
              data.frame(chr = seqnames(from),
                         start = rep(1, length(from)),
                         end = seqlengths(from)),
              ignore.strand = TRUE) %>%
              `names<-`(NULL)
      })

## Allow BSgenome to be converted to Views via GRanges.
setAs("BSgenome", "Views",
      function(from) as(from, "GRanges") %>% BSgenomeViews(from, .))

## Export the above coercion methods as described in Bioconductor S4 Objects lab
## exercise:
## https://www.bioconductor.org/help/course-materials/2011/AdvancedRFeb2011Seattle/ImplementingS4Objects-lab.pdf
## Here we use ROxygen tags to manage the NAMESPACE file for us instead of hand
## editing the file as suggested by the aged exercise.
#' @importFrom methods coerce
#' @exportMethod coerce
methods::coerce

#' Subset to standard DNA bases
#'
#' Functions for subsetting to standard DNA for a BSgenome (or BSgenomeViews).
#'
#' @param bsgenome The \code{\link[BSgenome]{BSgenome}} DNA to subset.
#' @param views The \code{\link[BSgenome]{BSgenomeViews}} DNA to subset.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation,
#'     or a \code{\link[base]{list}} of
#'     \code{\link[BiocParallel]{BiocParallelParam}} instances, to be applied in
#'     sequence for nested calls to \code{BiocParallel} functions.
#' @return A \code{BSgenomeViews} object with non-standard DNA bases (i.e. bases
#'     not in \code{\link[Biostrings]{DNA_BASES}}) removed for
#'     \code{\link{stddna_from_genome}}, and a \code{GRanges} object for
#'     \code{\link{stddna_from_views}}.
#' @name stddna
NULL

#' @rdname stddna
stddna_from_genome <- function(bsgenome, BPPARAM = bpparam()) {
    ## Convert to Views here for convenience, instead of subsetting by
    ## chromosomes.  Views allow more vector operations and preserve metadata
    ## better than shuttling around chromosomes.
    bsgenome %>% as("Views") %>% stddna_from_views(BPPARAM = BPPARAM) %>%
        BSgenomeViews(subject = bsgenome)
}

#' @rdname stddna
stddna_from_views <- function(views, BPPARAM = bpparam()) {
    bases_nonstd <- setdiff(DNA_ALPHABET, DNA_BASES)
    ## FIXME maybe parallelize this for larger genomes.
    bplapply(bases_nonstd, gr_masked, views = views, BPPARAM = BPPARAM) %>%
        sapply(gaps) %>% List() %>% unlist() %>% reduce() %>% gaps()
}

#' Remove motifs from BSgenomeViews
#'
#' Filter out DNA motifs to create a \code{\link[GenomicRanges]{GRanges-class}}
#' instance.
#'
#' @inheritParams stddna
#' @param motif The motif to remove in the sequences.
#' @return A \code{\link[GenomicRanges]{GRanges-class}} object with \code{motif}
#'     sequences removed.
gr_masked <- function(views, motif = "N") {
    views %>% as("DNAStringSet") %>% sapply(mask, motif) %>%
        as("RangesList") %>% ir2gr(views)
}

#' Convert IRanges back to GRanges using metadata from Views.
#'
#' We lose the \code{\link[GenomeInfoDb]{seqnames}} and
#' \code{\link[GenomeInfoDb]{seqinfo}} metadata when transforming the underlying
#' \code{\link[Biostrings]{DNAString}} of a
#' \code{\link[BSgenome]{BSgenomeViews}} object into
#' \code{\link[IRanges]{IRanges}}.  \code{\link{ir2gr}} recovers that
#' information from the original \code{\link[BSgenome]{BSgenomeViews}} to up
#' convert \code{\link[IRanges]{IRanges}} to
#' \code{\link[GenomicRanges]{GRanges}}.
#'
#' @inheritParams stddna
#' @param ranges IRanges
#' @return GRanges
#' @examples
#' bsgenome <- BSgenome::getBSgenome("BSgenome.Ecoli.NCBI.20080805")
#' views <- BSgenomeViews(bsgenome, GRanges("AC_000091:11-20"))
#' views
#' remaining <- mask(views[[1]], "T")
#' ir <- as(remaining, "IRanges")
#' ir
#' ## Convert IRanges to GRanges
#' ir2gr(ir, views)
ir2gr <- function(ranges, views) {
    ranges %>% shift(start(views) - 1) %>%
         `names<-`(seqnames(views)) %>%
         GRanges(seqinfo = seqinfo(views))
}

#' Chop up BSgenomeViews into overlapping k-mers.
#'
#' @param views The \code{\link[BSgenome]{BSgenomeViews}} sequences to segment
#'     into overlapping k-mers.
#' @param kmer The size of overlapping segments.
#' @return A \code{\link[BSgenome]{BSgenomeViews}} object with
#'     \code{\link[IRanges]{width}} exactly equal to \code{kmer} size.
#'     Ranges smaller than the \code{kmer} are dropped.
#' @export
kmerize <- function(views, kmer = 36) {
    ## Validate inputs.
    if (! is.numeric(kmer) | kmer %% 1 != 0)
        stop(sQuote("kmer"), " must be an integer, not ", sQuote(kmer))
    if (kmer < 1)
        stop(sQuote("kmer"), " must be >= 1, not ", sQuote(kmer))
    gr <- slidingWindows(granges(views), kmer) %>% unlist()
    BSgenomeViews(subject(views), gr)
}

#' Find hits of kmers along genome.
#' 
#' Compare \code{\link[Biostrings]{XString-class}} object with the
#' dictionary of the kmerized genome, then vectorize across 
#' 
#' @param xstring The \code{\link[Biostrings]{XString-class}}
#' @param views The \code{\link[BSgenome]{BSgenomeViews}}
#' @param pdict The \code{\link[Biostrings]{PDict}} dictionary.
#' @param indices Only operate on these indices of the \code{views}
#' @return A \code{\link[IRanges]{CompressedNormalIRangesList-class}} with no repeating hits.
#' @return A \code{\link[GenomicRanges]{GRanges-class}} with no repeating hits.

#' @rdname range
range_hits <- function(xstring, pdict) {
    matches <- matchPDict(pdict, xstring)
    hits <- elementNROWS(matches) > 1
    matches %>% as("CompressedIRangesList") %>% .[hits] %>%
        unlist() %>% IRanges::reduce()
}
#' @rdname range
ranges_hits <- function(views, pdict, indices = NULL) {
    if (is.null(indices)) indices <- 1:length(views)
    views_sub <- views[indices]
    irl <- lapply(views_sub, range_hits, pdict) %>% List()
    ## matchPDict() converted the BSgenomeViews into DNAStringSets, and
    ## therefore the start offsets and seqinfo was lost.  We need to readd that
    ## information using ir2gr() before combining all the ranges.
    lapply(seq_along(irl), function(i) ir2gr(irl[i], views_sub[i])) %>%
        List() %>% unlist() %>% GenomicRanges::reduce()
}

#' Return GRanges of uniquely mapping hits.
#'
#' @param views The \code{\link[BSgenome]{BSgenomeViews}} DNA to mapped.
#' @param genome The \code{\link[BSgenome]{BSgenome}} DNA to search for hits.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation,
#'     or a \code{\link[base]{list}} of
#'     \code{\link[BiocParallel]{BiocParallelParam}} instances, to be applied in
#'     sequence for nested calls to \code{BiocParallel} functions.
#' @return The \code{\link[GenomicRanges]{GRanges-class}} of mappable DNA sequences.
align <- function(views, genome = NULL, ...) {
    if (is.null(genome))
        genome <- subject(views)@pkgname
    file_fasta <- tempfile("views-", fileext = ".fasta")
    file_sample <- tempfile("sample-", fileext = ".txt")
    on.exit(unlink(c(file_fasta, file_sample)))
    writeXStringSet(as(views, "XStringSet") %>% unique(), file_fasta)
    write.table(data.frame(FileName = file_fasta,
                           SampleName = "kmers"), file_sample, sep = "\t",
                row.names = FALSE, quote = FALSE)
    ## Align!
    proj <- qAlign(file_sample, genome, ...)
    alignments(proj)[[1]]$FileName %>% readGAlignments() %>%
        `seqinfo<-`(value = seqinfo(views)) %>% `strand<-`(value = "*") %>%
        as("GRanges")
}

## FIXME yes, this function will be cleaned up.
#' Return mappable regions.
#' 
#' @param genome The \code{\link[BSgenome]{BSgenome}} DNA to be subset.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation,
#'     or a \code{\link[base]{list}} of
#'     \code{\link[BiocParallel]{BiocParallelParam}} instances, to be applied in
#'     sequence for nested calls to \code{BiocParallel} functions.
#' @return The \code{\link[GenomicRanges]{GRanges-class}} of mappable DNA sequences.
#' @examples
#' 
#' \dontrun{
#' mappable(BSgenome.Hsapiens.UCSC.hg38,
#' }
#'  
#' @export
mappable <- function(genome, BPPARAM = bpparam()) {
    message(sprintf(
        "[%6.2f sec] to load the genome",
        system.time(
            bsgenome <-
                BSgenome::getBSgenome(genome))[3]))
    ## ## Sanity check - there should be some non-DNA bases
    ## alphabetFrequency(as(bsgenome, "Views"), baseOnly = TRUE)
    message(sprintf(
        "[%6.2f sec] to subset to standard DNA bases in the genome",
        system.time(views <- stddna_from_genome(bsgenome, BPPARAM))[3]))
    message(sprintf(
        "[%6.2f sec] to chop up the genome into 36-mers",
        system.time(kmers <- kmerize(views))[3]))
    ## ## Sanity check - there should be no non-standard bases,
    ## ## otherwise pdict creation will fail.
    ## alphabetFrequency(kmers, baseOnly = TRUE) %>% colSums()
    ## ## There isn't really much one can do to speed up the PDict creation.
    message(sprintf(
        "[%6.2f sec] to create the search dictionary",
        system.time(pdict <- PDict(as(kmers,
                                      "DNAStringSet")))[3]))
    ## Chop up the views into smaller pieces.
    message(sprintf(
        "[%6.2f sec] to search the genome",
        system.time(gr <- bplapply(seq_along(views),
                                   ranges_hits,
                                   views = views,
                                   pdict = pdict,
                                   BPPARAM = BPPARAM) %>%
                        List() %>%
                        unlist() %>%
                        GenomicRanges::reduce())[3]))
    gr
}
