#' @importFrom BSgenome BSgenome BSgenomeViews subject elementNROWS
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom Biostrings DNAStringSet DNA_ALPHABET DNA_BASES PDict
#'     writeXStringSet matchPDict mask
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb seqinfo seqinfo<- seqlengths
#' @importFrom GenomicRanges GRanges granges end gaps reduce seqnames shift
#'     slidingWindows start strand<- width
#' @importFrom IRanges IRanges
#' @importFrom QuasR alignments qAlign
#' @importFrom S4Vectors List runLength<-
#' @importFrom magrittr %>%
#' @importFrom methods as
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
    ## matchPDict converted the BSgenomeViews into DNAStringSets, and
    ## therefore the start offsets and seqinfo was lost.  We need to readd that
    ## information using ir2gr before combining all the ranges.
    lapply(seq_along(irl), function(i) ir2gr(irl[i], views_sub[i])) %>%
        List() %>% unlist() %>% GenomicRanges::reduce()
}

#' Return GRanges of uniquely mapping hits.
#'
#' @param views The \code{\link[BSgenome]{BSgenomeViews}} DNA to mapped.
#' @param genome The \code{\link[BSgenome]{BSgenome}} DNA to search for hits.
#' @param ... Extra arguments passed on to \code{\link[QuasR]{qAlign}}.
#' @return The \code{\link[GenomicRanges]{GRanges-class}} of uniquely mapping DNA sequences.
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
            message(sprintf(" %.2f secs", time[3]))
        }
    }
    invisible(time[3])
}

#' Return mappable regions of a genome.
#'
#' @inheritParams kmerize
#' @param genome The \code{\link[BSgenome]{BSgenome}} DNA to be subset.
#' @param verbose Describe the calculation steps with time.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation,
#'     or a \code{\link[base]{list}} of
#'     \code{\link[BiocParallel]{BiocParallelParam}} instances, to be applied in
#'     sequence for nested calls to \code{BiocParallel} functions.
#' @return The \code{\link[GenomicRanges]{GRanges-class}} of mappable DNA sequences.
#' @examples
#' 
#' \dontrun{
#' mappable("BSgenome.Hsapiens.UCSC.hg38")
#' }
#'  
#' @export
mappable <- function(genome, kmer = 36, BPPARAM = bpparam(), verbose = TRUE) {
    timeit("Loading the genome", {
        bsgenome <- BSgenome::getBSgenome(genome)
    })
    timeit("Subsetting to standard DNA bases", {
        views <- stddna_from_genome(bsgenome, BPPARAM)
    })
    timeit("Chop up the genome into k-mers", {
        kmers <- kmerize(views, kmer = kmer)
    })
    timeit("Search the genome for unique k-mer alignments", {
        gr <- align(views, genome)
    })
    gr %>% GenomicRanges::reduce()
}
