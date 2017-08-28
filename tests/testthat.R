## Unset R_TESTS as workaround for QuasR starting another R process.
## https://stackoverflow.com/a/27994299 # nolint
Sys.setenv("R_TESTS" = "")
## Workaround for BiocParallel multicore failing with Travis CI, issue 7052:
## https://github.com/travis-ci/travis-ci/issues/7052 # nolint
BiocParallel:: register(BiocParallel::SerialParam())

library(testthat)
library(kmap)

test_check("kmap")
