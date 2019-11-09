library(testthat)
library(bayesplot)

Sys.unsetenv("R_TESTS")
test_check("bayesplot")
