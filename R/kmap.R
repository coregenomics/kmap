#' kmap: Generate mappable regions of the genome for a given K-mer length
#'
#' kmap mainly provides the \code{\link{mappable}} function.  It also provides
#' the utility function \code{\link{timeit}}.
#'
#' \code{\link{mappable}} uses Bowtie to determine regions with unique mapping
#' reads of a given length k; in other words, a given k-mer length.
#'
#' \code{\link{timeit}} prints messages about how long a code block takes to
#' execute.  It was originally created to provide informative, but suppressable,
#' messages for the code blocks in \code{mappable} which can take a long time to
#' run, and is generally useful for any long running function.
#'
#' @docType package
#' @name kmap
NULL
