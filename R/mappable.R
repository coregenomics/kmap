#' @importFrom BSgenome BSgenome BSgenomeViews subject
#' @importFrom BiocFileCache BiocFileCache bfccount bfcnew bfcquery bfcrpath
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom Biostrings DNAStringSet DNA_ALPHABET DNA_BASES writeXStringSet
#'     mask
#' @importFrom futile.logger flog.info
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb seqinfo seqinfo<- seqlengths
#' @importFrom GenomicRanges GRanges granges end gaps reduce seqnames shift
#'     slidingWindows start strand<-
#' @importFrom IRanges IRanges findOverlaps
#' @importFrom QuasR alignments qAlign
#' @importFrom S4Vectors List Rle elementNROWS from mendoapply runLength
#'     runLength<- runValue to
#' @importFrom methods as is
#' @importFrom plyr .
#' @importFrom utils write.table
NULL

## Reexport the `setAs` coercion methods below as described in Bioconductor S4
## Objects lab exercise:
## https://www.bioconductor.org/help/course-materials/2011/AdvancedRFeb2011Seattle/ImplementingS4Objects-lab.pdf # nolint
## Here we use ROxygen tags to manage the NAMESPACE file for us instead of hand
## editing the file as suggested by the above aged exercise.  The `setAs`
## functions themselves are intentionally not documented so as to not clobber
## the base R `as` help.  This creates the following check warning:
##
##   Undocumented S4 methods:
##     generic 'coerce' and siglist 'BSgenome,BSgenomeViews'
##     generic 'coerce' and siglist 'BSgenome,GRanges'
##     generic 'coerce' and siglist 'BSgenome,Views'
##   All user level objects in a package (including S4 classes and methods)
##
## Bioconductor core packages document `as` alongside the class descriptions or
## constructors, so these coercion functions should be upstreamed.
#' @importFrom methods coerce
#' @exportMethod coerce
methods::coerce
## Allow BSgenome to be converted to GRanges and Views object.
setAs("BSgenome", "GRanges",
      function(from) {
          gr <- GenomicRanges::makeGRangesFromDataFrame(
              data.frame(chr = seqnames(from),
                         start = rep(1, length(from)),
                         end = seqlengths(from)),
              ignore.strand = TRUE)
          names(gr) <- NULL
          gr
      })
## Allow BSgenome to be converted to BSgenomeViews via GRanges.
setAs("BSgenome", "BSgenomeViews",
      function(from) BSgenomeViews(from, as(from, "GRanges")))
## Alias Views to BSgenomeViews.
setAs("BSgenome", "Views",
      function(from) as(from, "BSgenomeViews"))

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
    if (! is(bsgenome, "BSgenomeViews"))
        stop(sQuote(bsgenome), "should be a BSgenomeViews object.")
    views <- stddna_from_views(bsgenome, BPPARAM = BPPARAM)
    BSgenomeViews(views, subject = subject(bsgenome))
}

#' @rdname stddna
stddna_from_views <- function(views, BPPARAM = bpparam()) {
    bases_nonstd <- setdiff(DNA_ALPHABET, DNA_BASES)
    ## FIXME maybe parallelize this for larger genomes.
    grl <- bplapply(bases_nonstd, gr_masked, views = views, BPPARAM = BPPARAM)
    grl <- List(sapply(grl, gaps))
    gr <- reduce(unlist(grl))
    gaps(gr)
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
    dna <- as(views, "DNAStringSet")
    rl <- as(sapply(dna, mask, motif), "IRangesList")
    rl <- mendoapply(shift, rl, start(views) - 1)
    seqnames <- expand_rle(seqnames(views), Rle(elementNROWS(rl)))
    strand <- "*"
    if (length(rl) == 0) {
        strand <- NULL
        seqnames <- NULL
    }
    GRanges(ranges = unlist(rl),
            strand = strand,
            seqnames = seqnames,
            seqinfo = seqinfo(views))
}

#' Insert into Rle runLength.
#'
#' @param rle Rle object to expand
#' @param size Rle with runValues of 1 for no insertion and runValues of > 1
#'     corresponding to insertions.  The runLength encodes the positions of
#'     insertions in \code{rle}
#' @return Rle with insertions
#' @examples
#' rle <- S4Vectors::Rle(c("a", "a", "a", "b", "c", "c"))
#' rle
#' size <- S4Vectors::Rle(c(1, 2, 1, 2, 2, 1))
#' kmap:::expand_rle(rle, size)
expand_rle <- function(rle, size) {
    ## Putting in the effort to modify the Rle is necessary because we're
    ## dealing with large memory objects.
    stopifnot(is(rle, "Rle"))
    stopifnot(is(size, "Rle"))
    stopifnot(length(rle) == length(size))
    ## Return original input if no expansion needed.
    if (all(size == 1)) {
        return(rle)
    }
    expand_by <- runValue(size) - 1
    to_expand <- expand_by > 0
    rle_ranges <- IRanges(start(rle), end(rle))
    ## Deduplicate run values by also checking end, not just start.
    start <- unlist(mapply(seq,
                           start(size)[to_expand],
                           end(size)[to_expand]))
    insertions <- IRanges(start = start, width = 1)
    ol <- findOverlaps(rle_ranges, insertions)
    expansion <- rep(0, length(rle_ranges))
    expansion[from(ol)] <- (size[start] - 1)[to(ol)]
    runLength(rle) <- runLength(rle) + expansion
    rle
}

#' Chop up BSgenomeViews into overlapping k-mers.
#'
#' @inheritParams stddna
#' @param views The \code{\link[BSgenome]{BSgenomeViews}} sequences to segment
#'     into overlapping k-mers.
#' @param kmer The size of overlapping segments.
#' @return A \code{\link[BSgenome]{BSgenomeViews}} object with
#'     \code{\link[IRanges]{width}} exactly equal to \code{kmer} size.
#'     Ranges smaller than the \code{kmer} are dropped.
kmerize <- function(views, kmer = 36) {
    ## Validate inputs.
    if (! is.numeric(kmer) | kmer %% 1 != 0)
        stop(sQuote("kmer"), " must be an integer, not ", sQuote(kmer))
    if (kmer < 1)
        stop(sQuote("kmer"), " must be >= 1, not ", sQuote(kmer))
    grl <- slidingWindows(granges(views), kmer)
    BSgenomeViews(subject(views), unlist(grl))
}

#' Return GRanges of uniquely mapping hits.
#'
#' @inheritParams stddna
#' @param views The \code{\link[BSgenome]{BSgenomeViews}} DNA to be mapped.
#' @param genome The \code{\link[BSgenome]{BSgenome}} DNA to search for hits.
#' @param ... Extra arguments passed on to \code{\link[QuasR]{qAlign}}.
#' @return The \code{\link[GenomicRanges]{GRanges-class}} of uniquely mapping DNA sequences.
align <- function(views, genome = NULL, BPPARAM = bpparam(), ...) {
    if (is.null(genome)) {
        ## nocov start
        genome <- subject(views)@pkgname
        ## nocov end
    } else if (is(genome, "BSgenomeViews")) {
        file_genome <- tempfile("genome-", fileext = ".fasta")
        on.exit(unlink(file_genome))
        xstringset <- as(genome, "XStringSet")
        names(xstringset) <- seqnames(genome)
        writeXStringSet(xstringset, file_genome)
        genome <- file_genome
    }
    ## QuasR does not yet support BiocParallel.
    cluster <- parallel::makeCluster(BPPARAM$workers)
    file_fasta <- tempfile("views-", fileext = ".fasta")
    file_sample <- tempfile("sample-", fileext = ".txt")
    on.exit({
        unlink(c(file_fasta, file_sample))
        parallel::stopCluster(cluster)
    }
    , add = TRUE)
    xstringset <- unique(as(views, "XStringSet"))
    writeXStringSet(xstringset, file_fasta)
    write.table(data.frame(FileName = file_fasta,
                           SampleName = "kmers"), file_sample, sep = "\t",
                row.names = FALSE, quote = FALSE)
    ## Align!
    proj <- qAlign(file_sample, genome, clObj = cluster, ...)
    gal <- readGAlignments(alignments(proj)[[1]]$FileName)
    seqinfo(gal) <- seqinfo(views)
    strand(gal) <- "*"
    as(gal, "GRanges")
}

#' Internal helper functions to manage mappable GRanges in BiocFileCache.
#'
#' Reference name for mappable genome of a given k-mer size.
#'
#' @name mappable_cache
#' @param bsgenome Apply the sequence information from this
#'     \code{\link[BSgenome]{BSgenome}} or \code{\link[BSgenome]{BSgenomeViews}}
#'     object to loaded \code{\link[GenomicRanges]{GRanges}}.
#' @inheritParams mappable
#' @return \code{\link{mappable_cache_name}} returns a character vector
#'     reference name for a given bsgenome and kmer size.
mappable_cache_name <- function(bsgenome, kmer = 36) {
    suffix <- NULL
    if (is(bsgenome, "BSgenome")) {
        pkgname <- bsgenome@pkgname
    } else if (is(bsgenome, "BSgenomeViews")) {
        pkgname <- subject(bsgenome)@pkgname
        ## If views subset the genome, add checksum to identifier.
        gr <- granges(bsgenome)
        gr_full <- as(subject(bsgenome), "GRanges")
        if (! identical(gr, gr_full)) {
            suffix <- digest::sha1(as.character(gr))
        }
    }
    words <- c("kmap", kmer, pkgname)
    if (! is.null(suffix))
        words <- c(words, suffix)
    paste(words, collapse = "_")
}

#' Initialize BiocFileCache assuming NULL to be default.
#'
#' @rdname mappable_cache
#' @return \code{\link{mappable_cache_bfc}} returns a
#'     \code{\link[BiocFileCache]{BiocFileCache}} object.
mappable_cache_bfc <- function(cache_path = NULL) {
    args <- list()
    if (! is.null(cache_path))
        args$cache <- cache_path
    do.call(BiocFileCache, args)
}

#' Path to cached GRanges.
#'
#' @rdname mappable_cache
#' @param name Reference name generated by \code{\link{mappable_cache_name}} to
#'     access \code{\link[BiocFileCache]{BiocFileCache}}.
#' @return \code{\link{mappable_cache_path}} returns a character vector of the
#'     path to the cached file path if the file exists, otherwise NULL.
mappable_cache_path <- function(name, cache_path = NULL) {
    name                                # Raise error if name was omitted.
    bfc <- mappable_cache_bfc(cache_path)
    query <- bfcquery(bfc, name)
    path <- NULL
    if (bfccount(query) == 1) {
        path <- bfcrpath(bfc, name)
    } else if (bfccount(query) > 1) {
        ## Name collision.
        stop("Found more than 1 cached file for", name)
    }
    path
}

#' Save GRanges to cache.
#'
#' @rdname mappable_cache
#' @param granges \code{\link[GenomicRanges]{GRanges}} to save to cache.
#' @param force Save to cache even if an object already exists.  Only used to
#'     simulate a name collision for unit tests and not intended for use
#'     otherwise; using this will break your ability to run
#'     \code{\link{mappable_cache_load}}!
#' @return \code{\link{mappable_cache_save}} returns the path to the saved cache
#'     file.
mappable_cache_save <- function(granges, name, cache_path = NULL,
                                force = FALSE) {
    path <- mappable_cache_path(name, cache_path)
    if (! force && ! is.null(path))
        return(path)
    bfc <- mappable_cache_bfc(cache_path)
    path <- suppressWarnings(bfcnew(bfc, name, ext = ".gff3"))
    rtracklayer::export.gff3(granges, path)
    path
}

#' Load GRanges from cache.
#'
#' @rdname mappable_cache
#' @return \code{\link{mappable_cache_load}} returns the
#'     \code{\link[GenomicRanges]{GRanges}} of mappable areas of the genome from
#'     the cache.
mappable_cache_load <- function(name, bsgenome, cache_path = NULL) {
    path <- mappable_cache_path(name, cache_path)
    if (is.null(path))
        return(NULL)
    gr <- rtracklayer::import(path)
    ## GFF3 does not serialize seqinfo.
    seqinfo(gr) <- seqinfo(bsgenome)
    gr
}

#' Return mappable regions of a genome for a given k-mer size.
#'
#' Genomes contain repetitive DNA, causing even long segments of DNA to map in
#' multiple places.  For a fixed segment size \code{kmer}, \code{link{mappable}}
#' checks whether DNA can map in multiple places in the genome and returns only
#' regions which uniquely map.
#'
#' @inheritParams kmerize
#' @inheritParams stddna
#' @param genome Single character vector genome identifier.  Short hand
#'     identifiers like "hg38" can be used in interactive sessions but can fail
#'     in scripts due to intended behavior of
#'     \code{\link[BSgenome]{getBSgenome}} preferring unambiguous full
#'     identifiers like "BSgenome.Hsapiens.UCSC.hg38".
#' @param kmer The exact size of DNA segments to check for unique mapping in the
#'     genome.
#' @param verbose Print messages of the calculation steps and their elapsed
#'     time.  Unfortunately, during first-time calculation, messages from
#'     \code{\link[QuasR]{qAlign}} are always shown because \code{qAlign} runs
#'     as a subprocess.
#' @param cache_path Path to initialize
#'     \code{\link[BiocFileCache]{BiocFileCache}} object at a non-default path
#'     storage location.
#' @return The \code{\link[GenomicRanges]{GRanges-class}} of mappable DNA
#'     sequences.
#' @examples
#' \dontrun{
#' mappable("BSgenome.Scerevisiae.UCSC.sacCer2")
#' }
#'  
#' @export
mappable <- function(genome, kmer = 36, BPPARAM = bpparam(), verbose = TRUE,
                     cache_path = NULL) {
    if (is(genome, "BSgenomeViews")) {
        ## Use potentially subset genome.
        bsgenome <- genome
    } else if (is.character(genome)) {
        ## Validate genome string.
        bsgenome <- BSgenome::getBSgenome(genome)
    } else {
        stop(sQuote(genome),
             "must be a genome identifier or BSgenomeViews object")
    }
    ## Return the cached mappable object if it already exists.
    name <- mappable_cache_name(bsgenome, kmer)
    path <- mappable_cache_path(name, cache_path)
    if (! is.null(path)) {
        flog.info(sprintf("Reading the cached mappable genome"))
        gr <- mappable_cache_load(name, bsgenome, cache_path)
        return(gr)
    }
    ## No cached object available, so calculate mappable GRanges.
    flog.info("Removing non-standard DNA bases")
    views <- stddna_from_genome(as(bsgenome, "BSgenomeViews"), BPPARAM)
    flog.info(sprintf("Chopping into %d-mers", kmer))
    views <- kmerize(views, kmer = kmer)
    flog.info(sprintf("Search the genome for unique %d-mer alignments", kmer))
    gr <- reduce(align(views, genome, BPPARAM))
    ## Suppress warnings from BiocFileCache.
    flog.info("Saving result to cache")
    mappable_cache_save(gr, name, cache_path)
    gr
}
