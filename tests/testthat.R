## Unset R_TESTS as workaround for QuasR starting another R process.
## https://stackoverflow.com/a/27994299 # nolint
Sys.setenv("R_TESTS" = "")

library(testthat)
library(kmap)

test_check("kmap")
