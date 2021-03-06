% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mappable.R
\name{align}
\alias{align}
\title{Return GRanges of uniquely mapping hits.}
\usage{
align(views, genome = NULL, BPPARAM = bpparam(), ...)
}
\arguments{
\item{views}{The \code{\link[BSgenome]{BSgenomeViews}} DNA to be mapped.}

\item{genome}{The \code{\link[BSgenome]{BSgenome}} DNA to search for hits.}

\item{BPPARAM}{An optional \code{\link[BiocParallel]{BiocParallelParam}}
instance determining the parallel back-end to be used during evaluation,
or a \code{\link[base]{list}} of
\code{\link[BiocParallel]{BiocParallelParam}} instances, to be applied in
sequence for nested calls to \code{BiocParallel} functions.}

\item{...}{Extra arguments passed on to \code{\link[QuasR]{qAlign}}.}
}
\value{
The \code{\link[GenomicRanges]{GRanges-class}} of uniquely mapping DNA sequences.
}
\description{
Return GRanges of uniquely mapping hits.
}
