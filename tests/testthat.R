## R CMD Check workarounds:
##
## Unset R_TESTS for QuasR starting another R process:
## https://stackoverflow.com/a/27994299 # nolint
Sys.setenv("R_TESTS" = "")
## Travis CI workarounds:
##
## BiocParallel multicore fails:
## https://github.com/travis-ci/travis-ci/issues/7052 # nolint
BiocParallel:: register(BiocParallel::SerialParam())
## Don't treat warnings as errors, because 2 warnings are non-trivial to fix.
options(warn = 1)

library(testthat)
library(kmap)

test_check("kmap")
