## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
Sys.setenv(USE_CXX14 = "1")

## ------------------------------------------------------------------------
optim(rnorm(3), fn = f, gr = g, a = c(1, 2, 3), method = "BFGS")$par  # Rcpp exported f and g

