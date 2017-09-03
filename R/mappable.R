#' @importFrom BSgenome BSgenome BSgenomeViews subject
#' @importFrom BiocFileCache BiocFileCache bfccount bfcnew bfcquery bfcrpath
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom Biostrings DNAStringSet DNA_ALPHABET DNA_BASES writeXStringSet
#'     mask
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb seqinfo seqinfo<- seqlengths
#' @importFrom GenomicRanges GRanges granges end gaps reduce seqnames shift
#'     slidingWindows start strand<-
#' @importFrom IRanges IRanges
#' @importFrom QuasR alignments qAlign
#' @importFrom S4Vectors List
#' @importFrom magrittr %>%
#' @importFrom methods as is
#' @importFrom plyr .
#' @importFrom utils write.table
NULL

## Reexport to use in the the unittest suite, etc.
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

## Export the `setAs` coercion methods below as described in Bioconductor S4 Objects
## lab exercise:
## https://www.bioconductor.org/help/course-materials/2011/AdvancedRFeb2011Seattle/ImplementingS4Objects-lab.pdf # nolint
## Here we use ROxygen tags to manage the NAMESPACE file for us instead of hand
## editing the file as suggested by the above aged exercise.  The `setAs`
## functions themselves are intentionally not documented so as to not clobber
## the base R `as` help.  This creates the following check warning:
##
##   Undocumented S4 methods:
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
    bsgenome %>% stddna_from_views(BPPARAM = BPPARAM) %>%
        BSgenomeViews(subject = subject(bsgenome))
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
#' ## This example finds GRanges not containing the DNA base "T"
#' ## in a section of the Ecoli genome.
#' bsgenome <- BSgenome::getBSgenome("BSgenome.Ecoli.NCBI.20080805")
#' views <- BSgenomeViews(bsgenome, GRanges("AC_000091:11-20"))
#' views
#' remaining <- mask(views[[1]], "T")
#' ir <- as(remaining, "IRanges")
#' ir
#' ## ir now has more rows than views after being split by `mask`.
#' ## Duplicate the rows to match.
#' views <- rep(views, length.out = length(ir))
#' ## Convert IRanges to GRanges
#' kmap:::ir2gr(ir, views)
ir2gr <- function(ranges, views) {
    if (NROW(ranges) != NROW(views))
        stop("Length of IRanges must match Views")
    if (class(ranges) != "RangesList")
        ranges <- as(ranges, "RangesList")
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
kmerize <- function(views, kmer = 36) {
    ## Validate inputs.
    if (! is.numeric(kmer) | kmer %% 1 != 0)
        stop(sQuote("kmer"), " must be an integer, not ", sQuote(kmer))
    if (kmer < 1)
        stop(sQuote("kmer"), " must be >= 1, not ", sQuote(kmer))
    gr <- slidingWindows(granges(views), kmer) %>% unlist()
    BSgenomeViews(subject(views), gr)
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
        writeXStringSet(as(genome, "XStringSet") %>%
                        `names<-`(seqnames(genome)),
                        file_genome)
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
    writeXStringSet(as(views, "XStringSet") %>% unique(), file_fasta)
    write.table(data.frame(FileName = file_fasta,
                           SampleName = "kmers"), file_sample, sep = "\t",
                row.names = FALSE, quote = FALSE)
    ## Align!
    proj <- qAlign(file_sample, genome, clObj = cluster, ...)
    alignments(proj)[[1]]$FileName %>% readGAlignments() %>%
        `seqinfo<-`(value = seqinfo(views)) %>% `strand<-`(value = "*") %>%
        as("GRanges")
}

#' Run code with message and timing information.
#'
#' Only print messages if the magic variable \code{verbose} is set in the
#' environment from where \code{timeit} was called.
#'
#' @param msg Message describing the code block being run.
#' @param code Block of code to evaluate and time.  Enclose in curly braces for
#'     multiple lines.
#' @return Invisibly returns evaluation time taken.
#' @examples
#' ## Use timeit directly inside a script
#' verbose <- TRUE
#' timeit("My long calculation", Sys.sleep(0.5))
#' ## Silence the same calculation
#' verbose <- FALSE
#' sec <- timeit("My long calculation", Sys.sleep(0.5))
#' sec
#'
#' ## Use timeit inside a function
#' my_long_calculation <- function(verbose = TRUE) {
#'     timeit("Step 1", Sys.sleep(0.1))
#'     timeit("Step 2", {
#'         x <- 10
#'         Sys.sleep(0.2)
#'     })
#'     x
#' }
#' my_long_calculation()
#' my_long_calculation(verbose = FALSE)
#' 
#' @export
timeit <- function(msg, code) {
    env <- parent.frame()
    verbose <- FALSE
    ## Check if verbose was defined in the parent environment.
    if (exists("verbose", envir = env))
        verbose <- env$verbose
    if (verbose) message(msg, " ...", appendLF = FALSE)
    time <- system.time(eval(code, env))
    if (verbose) {
        ## Use the lubridate package for automatic rounding to minutes, hours,
        ## etc for long calculations.
        if (suppressPackageStartupMessages(
            requireNamespace("lubridate", quietly = TRUE))) {
            diff <- lubridate::make_difftime(time[3])
            message(sprintf(" %.2f %s", diff, attr(diff, "units")))
        } else {
            ## nocov start
            message(sprintf(" %.2f secs", time[3]))
            ## nocov end
        }
    }
    invisible(time[3])
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
            suffix <- gr %>% as.character() %>% digest::sha1()
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
        timeit(sprintf("Reading the cached mappable genome"), {
            gr <- mappable_cache_load(name, bsgenome, cache_path)
        })
        return(gr)
    }
    ## No cached object available, so calculate mappable GRanges.
    timeit(sprintf("Removing non-standard DNA bases and chopping into %d-mers",
                   kmer),
           views <- stddna_from_genome(bsgenome, BPPARAM) %>%
               kmerize(kmer = kmer))
    timeit(sprintf("Search the genome for unique %d-mer alignments", kmer),
           gr <- align(views, genome, BPPARAM) %>% reduce())
    ## Suppress warnings from BiocFileCache.
    timeit("Saving result to cache",
           mappable_cache_save(gr, name, cache_path))
    gr
}
