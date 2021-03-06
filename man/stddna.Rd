% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mappable.R
\name{stddna}
\alias{stddna}
\alias{stddna_from_genome}
\alias{stddna_from_views}
\title{Subset to standard DNA bases}
\usage{
stddna_from_genome(bsgenome, BPPARAM = bpparam())

stddna_from_views(views, BPPARAM = bpparam())
}
\arguments{
\item{bsgenome}{The \code{\link[BSgenome]{BSgenome}} DNA to subset.}

\item{BPPARAM}{An optional \code{\link[BiocParallel]{BiocParallelParam}}
instance determining the parallel back-end to be used during evaluation,
or a \code{\link[base]{list}} of
\code{\link[BiocParallel]{BiocParallelParam}} instances, to be applied in
sequence for nested calls to \code{BiocParallel} functions.}

\item{views}{The \code{\link[BSgenome]{BSgenomeViews}} DNA to subset.}
}
\value{
A \code{BSgenomeViews} object with non-standard DNA bases (i.e. bases
    not in \code{\link[Biostrings]{DNA_BASES}}) removed for
    \code{\link{stddna_from_genome}}, and a \code{GRanges} object for
    \code{\link{stddna_from_views}}.
}
\description{
Functions for subsetting to standard DNA for a BSgenome (or BSgenomeViews).
}
