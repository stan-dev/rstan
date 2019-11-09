params <-
list(EVAL = TRUE)

## ---- SETTINGS-knitr, include=FALSE--------------------------------------
stopifnot(require(knitr))
opts_chunk$set(
  comment=NA,
  eval = if (isTRUE(exists("params"))) params$EVAL else FALSE,
  dev = "png",
  dpi = 150,
  fig.asp = 0.618,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
)

## ---- setup, message=FALSE-----------------------------------------------
library("rstanarm")
library("bayesplot")
library("loo")

## ------------------------------------------------------------------------
# the 'roaches' data frame is included with the rstanarm package
data(roaches)
str(roaches)

# rescale to units of hundreds of roaches
roaches$roach1 <- roaches$roach1 / 100

## ---- count-roaches-mcmc, results="hide"---------------------------------
fit1 <-
  stan_glm(
    formula = y ~ roach1 + treatment + senior,
    offset = log(exposure2),
    data = roaches,
    family = poisson(link = "log"),
    prior = normal(0, 2.5, autoscale = TRUE),
    prior_intercept = normal(0, 5, autoscale = TRUE),
    seed = 12345
  )

## ------------------------------------------------------------------------
loo1 <- loo(fit1, save_psis = TRUE)

## ------------------------------------------------------------------------
print(loo1)

## ---- out.width = "70%"--------------------------------------------------
plot(loo1)

## ------------------------------------------------------------------------
yrep <- posterior_predict(fit1)

# requires bayesplot version >= 1.5.0
ppc_loo_pit_overlay(
  y = roaches$y,
  yrep = yrep,
  lw = weights(loo1$psis_object)
)

## ---- count-roaches-negbin, results="hide"-------------------------------
fit2 <- update(fit1, family = neg_binomial_2)

## ------------------------------------------------------------------------
loo2 <- loo(fit2, save_psis = TRUE, cores = 2)
print(loo2)

## ------------------------------------------------------------------------
plot(loo2, label_points = TRUE)

## ------------------------------------------------------------------------
if (any(pareto_k_values(loo2) > 0.7)) {
  loo2 <- loo(fit2, save_psis = TRUE, k_threshold = 0.7)
}
print(loo2)

## ------------------------------------------------------------------------
yrep <- posterior_predict(fit2)
ppc_loo_pit_overlay(roaches$y, yrep, lw = weights(loo2$psis_object))

## ---- count-roaches-loo--------------------------------------------------
compare_models(loo1, loo2)  # use loo::compare(loo1, loo2) if not using rstanarm

