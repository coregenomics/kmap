# kmap

[![Build Status](https://travis-ci.org/coregenomics/kmap.svg?branch=master)](https://travis-ci.org/coregenomics/kmap)
[![codecov.io](https://codecov.io/gh/coregenomics/kmap/branch/master/graphs/badge.svg)](https://codecov.io/gh/coregenomics/kmap)

Generate mappable regions of the genome for a given K-mer length.

Knowing where reads can uniquely map in the genome
is useful for nascent RNA assays,
both in statistical calculations and to make predictions.
kmap takes two inputs: a genome and the read length, K.
The output is automatically cached on disk.

``` R
library(kmap)

gr <- mappable("hg38", kmer = 36)
```

## Usage

Install the latest version of kmer using:

``` sh
devtools::install_github("coregenomics/kmap", repos = BiocInstaller::biocinstallRepos())
```

If the above command fails, install Bioconductor and the `devtools` package.

``` R
source("https://bioconductor.org/biocLite.R")
install.packages("devtools")
```
