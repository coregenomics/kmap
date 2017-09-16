#' kmap: Generate mappable regions of the genome for a given K-mer length
#'
#' kmap provides the \code{\link{mappable}} function.
#'
#' \code{\link{mappable}} uses Bowtie to determine regions with unique mapping
#' reads of a given length k; in other words, a given k-mer length.
#'
#' Finding mappable regions can take a few minutes even for small genomes,
#' therefore \code{\link{mappable}} prints informative messages about processing
#' steps using the \code{\link[futile.logger]{futile.logger}} package.
#'
#' @docType package
#' @name kmap
NULL
