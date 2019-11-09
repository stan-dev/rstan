params <-
list(EVAL = TRUE)

## ---- SETTINGS-knitr, include=FALSE--------------------------------------
stopifnot(require(knitr))
opts_chunk$set(
  comment=NA,
  eval = params$EVAL
)

## ---- eval=FALSE---------------------------------------------------------
#  library("rstantools")
#  rstan_package_skeleton(path = 'rstanlm')

## ---- echo=FALSE,warning=FALSE-------------------------------------------
library("rstantools")
td <- tempdir()
PATH <- file.path(td, "rstanlm")
rstan_package_skeleton(path = PATH, rstudio=FALSE, open=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  setwd("rstanlm")
#  list.files(all.files = TRUE)

## ---- echo=FALSE---------------------------------------------------------
list.files(PATH, all.files = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  file.show("DESCRIPTION")

## ---- echo=FALSE---------------------------------------------------------
DES <- readLines(file.path(PATH, "DESCRIPTION"))
cat(DES, sep = "\n")

## ---- eval=FALSE---------------------------------------------------------
#  file.show("Read-and-delete-me")

## ---- echo=FALSE---------------------------------------------------------
cat(readLines(file.path(PATH, "Read-and-delete-me")), sep = "\n")

## ---- eval=FALSE---------------------------------------------------------
#  file.remove('Read-and-delete-me')

## ---- echo=FALSE---------------------------------------------------------
file.remove(file.path(PATH, 'Read-and-delete-me'))

## ---- include=FALSE------------------------------------------------------
stan_prog <- "
data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real intercept;
  real beta;
  real<lower=0> sigma;
}
model {
  // ... priors, etc.
  
  y ~ normal(intercept + beta * x, sigma);
}
"
cat(stan_prog, file = file.path(PATH, "src", "stan_files", "lm.stan"))

## ------------------------------------------------------------------------
# Save this file as `R/lm_stan.R`

#' Bayesian linear regression with Stan
#'
#' @export
#' @param x Numeric vector of input values.
#' @param y Numberic vector of output values.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
lm_stan <- function(x, y, ...) {
  standata <- list(x = x, y = y, N = length(y))
  out <- rstan::sampling(stanmodels$lm, data = standata, ...)
  return(out)
}

## ---- include=FALSE------------------------------------------------------
Rcode <- "
#' Bayesian linear regression with Stan
#'
#' @export
#' @param x Numeric vector of input values.
#' @param y Numberic vector of output values.
#' @param ... Arguments passed to `rstan::sampling`.
#' @return An object of class `stanfit` returned by `rstan::sampling`
lm_stan <- function(x, y) {
  out <- rstan::sampling(stanmodels$lm, data=list(x=x, y=y, N=length(y)))
  return(out)
}
"
cat(Rcode, file = file.path(PATH, "R", "lm_stan.R"))

## ---- eval=FALSE---------------------------------------------------------
#  file.show(file.path("R", "rstanlm-package.R"))

## ---- echo=FALSE---------------------------------------------------------
cat(readLines(file.path(PATH, "R", "rstanlm-package.R")), sep = "\n")

## ---- eval=FALSE---------------------------------------------------------
#  roxygen2::roxygenise(clean=TRUE)

## ---- echo=FALSE---------------------------------------------------------
roxygen2::roxygenise(PATH, clean=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  devtools::install(local=FALSE)

## ----echo=FALSE, results="hide"------------------------------------------
devtools::load_all(PATH, recompile=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  library("rstanlm")

## ------------------------------------------------------------------------
fit <- lm_stan(y = rnorm(10), x = rnorm(10), iter = 500)
print(fit)

## ---- echo=FALSE---------------------------------------------------------
unlink(PATH, recursive = TRUE)

