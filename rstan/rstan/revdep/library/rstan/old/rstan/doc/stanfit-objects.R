## ----setup, include=FALSE------------------------------------------------
library(rstan)
knitr::opts_chunk$set(
  echo = TRUE,
  comment = NA,
  fig.align = "center",
  fig.height = 5,
  fig.width = 7
  )

## ---- example-model, eval=FALSE------------------------------------------
#  library(rstan)
#  fit <- stan_demo("eight_schools", refresh = 0)

## ---- fit, echo=FALSE, cache=FALSE, results="hide"-----------------------
J <- 8
y <- c(28,  8, -3,  7, -1,  1, 18, 12)
sigma <- c(15, 10, 16, 11,  9, 11, 10, 18)
fit <- stan(
  file= "schools.stan", 
  model_name = "eight_schools",
  data = c("y", "J", "sigma"), 
  refresh = 0
  )

## ---- stanfit-class------------------------------------------------------
class(fit)

## ---- extract-1----------------------------------------------------------
list_of_draws <- extract(fit)
print(names(list_of_draws))

## ---- extract-2----------------------------------------------------------
head(list_of_draws$mu)
head(list_of_draws$tau)
head(list_of_draws$theta)

## ---- as.matrix-1--------------------------------------------------------
matrix_of_draws <- as.matrix(fit)
print(colnames(matrix_of_draws))

df_of_draws <- as.data.frame(fit)
print(colnames(df_of_draws))

array_of_draws <- as.array(fit)
print(dimnames(array_of_draws))

## ---- as.matrix-2, results="hold"----------------------------------------
print(dim(matrix_of_draws))
print(dim(df_of_draws))
print(dim(array_of_draws))

## ---- as.matrix-3--------------------------------------------------------
mu_and_theta1 <- as.matrix(fit, pars = c("mu", "theta[1]"))
head(mu_and_theta1)

## ---- summary-1----------------------------------------------------------
fit_summary <- summary(fit)
print(names(fit_summary))

## ---- summary-2----------------------------------------------------------
print(fit_summary$summary)

## ---- summary-3----------------------------------------------------------
mu_tau_summary <- summary(fit, pars = c("mu", "tau"), probs = c(0.1, 0.9))$summary
print(mu_tau_summary)

## ---- summary-4----------------------------------------------------------
mu_tau_80pct <- mu_tau_summary[, c("10%", "90%")]
print(mu_tau_80pct)

## ---- get_sampler_params-1-----------------------------------------------
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
sampler_params_chain1 <- sampler_params[[1]]
colnames(sampler_params_chain1)

## ---- get_sampler_params-2-----------------------------------------------
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
print(mean_accept_stat_by_chain)

max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
print(max_treedepth_by_chain)

## ---- get_stan_code-1----------------------------------------------------
code <- get_stancode(fit)

## ---- get_stan_code-2----------------------------------------------------
print(code)

## ---- get_stan_code-3----------------------------------------------------
cat(code)

## ---- get_inits----------------------------------------------------------
inits <- get_inits(fit)
inits_chain1 <- inits[[1]]
print(inits_chain1)

## ---- get_seed-----------------------------------------------------------
print(get_seed(fit))

## ---- get_elapsed_time---------------------------------------------------
print(get_elapsed_time(fit))

