## ----setup, include=FALSE------------------------------------------------
library(rstan)
knitr::opts_chunk$set(
  echo = TRUE,
  comment = NA,
  fig.align = "center",
  fig.height = 5,
  fig.width = 7
  )

## ---- echo=FALSE, comment=""---------------------------------------------
cat(readLines("schools.stan"), sep = "\n")

## ---- lookup-------------------------------------------------------------
lookup("dnorm")
lookup(dwilcox)   # no corresponding Stan function

## ---- schools-data-------------------------------------------------------
schools_data <- list(
  J = 8,
  y = c(28,  8, -3,  7, -1,  1, 18, 12),
  sigma = c(15, 10, 16, 11,  9, 11, 10, 18)
)

## ---- callstan, results="hide"-------------------------------------------
library(rstan)
fit1 <- stan(
  file = "schools.stan",  # Stan program
  data = schools_data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  refresh = 0             # no progress shown
  )

## ---- print--------------------------------------------------------------
print(fit1, pars=c("theta", "mu", "tau", "lp__"), probs=c(.1,.5,.9))

## ---- stanfit-plot-------------------------------------------------------
plot(fit1)

## ---- stanfit-traceplot--------------------------------------------------
traceplot(fit1, pars = c("mu", "tau"), inc_warmup = TRUE, nrow = 2)

## ---- stanfit-print------------------------------------------------------
print(fit1, pars = c("mu", "tau"))

## ---- get_sampler_params-------------------------------------------------
# all chains combined
sampler_params <- get_sampler_params(fit1, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

# each chain separately
lapply(sampler_params, summary, digits = 2)

## ---- pairs-plot---------------------------------------------------------
pairs(fit1, pars = c("mu", "tau", "lp__"), las = 1)

## ---- expose_stan_functions----------------------------------------------
model_code <-
'
functions {
  real standard_normal_rng() {
    return normal_rng(0,1);
  }
}
model {}
'
expose_stan_functions(stanc(model_code = model_code))
standard_normal_rng()

## ---- optimizer, results="hide"------------------------------------------
ocode <- "
  data {
    int<lower=1> N;
    real y[N];
  } 
  parameters {
    real mu;
  } 
  model {
    target += normal_lpdf(y | mu, 1);
  } 
"

sm <- stan_model(model_code = ocode)
y2 <- rnorm(20)

## ------------------------------------------------------------------------
mean(y2)
optimizing(sm, data = list(y = y2, N = length(y2)), hessian = TRUE)

