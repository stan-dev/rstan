Sys.setenv(USE_CXX14 = 1)
library(rstan)

rstan_options(auto_write = TRUE)

## This is a demo for solving ODEs with an analytic Jacobian which is
## automatically generated.

set.seed(234235)

source("./integrate_jacobian_analytic/make_ode_model.R")

stan_file <- "oral_2cmt_jacobian_analytic.stan"

cpp_ode_dir <- normalizePath("./integrate_jacobian_analytic")

sm <- stan_ode_model(stan_file, cpp_ode_dir)

sm_par <- stan_ode_model(stan_file, cpp_ode_dir, enable_threads=TRUE)

## to generate the hpps and look at them, do (can be used with cmdstan
## as well)
sml <- stan_ode_model(stan_file, cpp_ode_dir, compile=FALSE)

## setup stan data

prior_theta_mean <- c(log(2)/c(1, 32, 8), 0.8) ## ka, CL, Q, V2
prior_theta_sd <- rep(0.25, 4)

ke <- prior_theta_mean[2]
k12 <- prior_theta_mean[3]
k21 <- prior_theta_mean[3] / prior_theta_mean[4]
D <- (k12 + k21 + ke)^2 - 4 * k21 * ke
D

rel_tol <- 1E-8
abs_tol <- 1E-5
max_steps <- 1E3
solver <- 0

T <- 2 * 24
J <- 6 * 5
data <- list(T=T, J=J, prior_theta_mean=prior_theta_mean, prior_theta_sd=prior_theta_sd,
             yobs=rep(1, J*T), rel_tol=rel_tol, abs_tol=abs_tol, max_steps=max_steps, solver=solver, analytic_jacobian=1)

## useful for cmdstan runs
stan_rdump(names(data), "oral_2cmt_jac_data.R", envir=list2env(data))

## simulate some truth which we recover later

## simulate using parallel and serial code

eta_v <- array(rnorm(J, 0, 1), dim=c(J,1))
truth <- list(theta=prior_theta_mean, eta=eta_v, sigma_y=0.1)

sim <- sampling(sm, data=data, warmup=0, iter=1, chains=1, seed=120,
                algorithm="Fixed", init=list(truth))

yhat_sim <- extract(sim, pars="yhat")[[1]][1,]
yrepl_sim <- extract(sim, pars="yrepl")[[1]][1,]

plot(1:T, yrepl_sim[1:T])
lines(1:T, yhat_sim[1:T])

plot(1:T, log(yrepl_sim[1:T]))
lines(1:T, log(yhat_sim[1:T]))

## now fit the simulated data set in parallel and serially, then compare

fake_data_rk45 <- modifyList(data, list(observed=1, yobs=yrepl_sim, solver=0))
fake_data_bdf <- modifyList(fake_data_rk45, list(solver=1))
fake_data_adams <- modifyList(fake_data_rk45, list(solver=2))

fake_data_rk45_ad <- modifyList(fake_data_rk45, list(analytic_jacobian=0))
fake_data_bdf_ad <- modifyList(fake_data_rk45, list(solver=1))
fake_data_adams_ad <- modifyList(fake_data_rk45, list(solver=2))


control <- list(adapt_delta=0.8, stepsize=1)

system.time(fake_run_rk45 <- sampling(sm, data=fake_data_rk45, warmup=250, iter=400, chains=1, init=list(truth), seed=12, refresh=10,
                                      control=control))

system.time(fake_run_rk45_ad <- sampling(sm, data=fake_data_rk45_ad, warmup=250, iter=400, chains=1, init=list(truth), seed=12, refresh=10,
                                         control=control))

## use as many cores as the machine has available
Sys.setenv(STAN_NUM_THREADS=-1)
Sys.getenv("STAN_NUM_THREADS")

system.time(fake_run_par_rk45 <- sampling(sm_par, data=fake_data_rk45, warmup=250, iter=400, chains=1, init=list(truth), seed=12, refresh=10,
                                      control=control))

system.time(fake_run_par_rk45_ad <- sampling(sm_par, data=fake_data_rk45_ad, warmup=250, iter=400, chains=1, init=list(truth), seed=12, refresh=10,
                                             control=control))

## now look at sampler parameters and sample quality
colMeans(get_sampler_params(fake_run_rk45, inc_warmup=FALSE)[[1]])
colMeans(get_sampler_params(fake_run_par_rk45, inc_warmup=FALSE)[[1]])

pars <- c("theta", "eta", "sigma_y")
print(fake_run_rk45, pars)
print(fake_run_par_rk45, pars)

truth[[1]]
truth[[2]]

post_serial <- as.matrix(fake_run_rk45, pars=pars)
post_par <- as.matrix(fake_run_par_rk45, pars=pars)

## compare estimates vs the truth
library(bayesplot)

true <- as.vector(unlist(truth))

bias_par <- sweep(post_par, 2, true)

mcmc_recover_intervals(bias_par, rep(0, length(true)))

bias_serial <- sweep(post_serial, 2, true)

mcmc_recover_intervals(bias_serial, rep(0, length(true)))
