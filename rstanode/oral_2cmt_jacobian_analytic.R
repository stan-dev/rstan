library(rstan)
rstan_options(auto_write = TRUE)

## This is a demo for solving ODEs with an analytic Jacobian which is
## automatically generated.

set.seed(234235)

source("./integrate_jacobian_analytic/make_ode_model.R")

stan_file <- "oral_2cmt_jacobian_analytic.stan"

cpp_ode_dir <- normalizePath("./integrate_jacobian_analytic")

sm <- stan_ode_model(stan_file, cpp_ode_dir)


## to generate the hpps and look at them, do
sml <- stan_ode_model(stan_file, cpp_ode_dir, compile=FALSE)

## setup stan data
state0 <- c(1000, 0, 0)
t0 <- -0.1
ts <- seq(0, 24*3, length=11)
theta <- log(2)/c(3, 12, 7, 10)
theta_sd <- rep(0.5, 4)
dose <- 10
xr <- c(1.)
xi <- c(1L)
rel_tol <- 1E-9
abs_tol <- 1E-6
max_steps <- 1E3
solver <- 1
## note: the initials are of order 1e2 while the parameters 1e-1
pbar <- c(rep(100, 3), rep(1, 4))

T <- 2 * 24
J <- 4 * 3
data_par <- list(T=T, J=J, theta=theta, theta_sd=theta_sd, dose=dose, parallel=1, observed=0, yobs=array(dim=0), rel_tol=rel_tol, abs_tol=abs_tol, max_steps=max_steps, solver=solver, pbar=pbar)
data_serial <- modifyList(data_par, list(parallel=0))

## useful for cmdstan runs
stan_rdump(names(data_par), "oral_2cmt_jac_data_par.R", envir=list2env(data_par))
stan_rdump(names(data_serial), "oral_2cmt_jac_data_serial.R", envir=list2env(data_serial))

Sys.getenv("OMP_NUM_THREADS")

## simulate some truth which we recover later

## simulate using parallel and serial code
truth <- list(theta_v=theta, dose0_v=dose + 1:J - 1, sigma_y=0.2)
sim <- sampling(sm, data=data_serial, warmup=0, iter=1, chains=1, seed=12,
                algorithm="Fixed", init=list(truth))

sim2 <- sampling(sm, data=data_par, warmup=0, iter=1, chains=1, seed=12,
                algorithm="Fixed", init=list(truth))


yhat_sim <- extract(sim, pars="yhat")[[1]][1,,]
yrepl_sim <- extract(sim, pars="yrepl")[[1]][1,]

yhat_sim2 <- extract(sim2, pars="yhat")[[1]][1,,]
yrepl_sim2 <- extract(sim2, pars="yrepl")[[1]][1,]

## ode values match? Only for bdf integrator and pbar==1 they match
## exactly
all(yhat_sim == yhat_sim2)
all(yrepl_sim == yrepl_sim2)

## does the gradient match?

truth_un <- unconstrain_pars(sim, truth)

log_prob(sim, truth_un)
grad_sim <- grad_log_prob(sim, truth_un)
grad_sim

log_prob(sim2, truth_un)
grad_sim2 <- grad_log_prob(sim2, truth_un)
grad_sim2

## note: we only get equivalence at the order of integration
## tolerances
all(grad_sim == grad_sim2)

grad_sim - grad_sim2

summary(yhat_sim)

library(ggplot2)
theme_set(theme_bw())

yhat_sim <- as.data.frame(yhat_sim)

yhat_sim <- transform(yhat_sim, time=1:T, id=rep(1:J, each=T))

names(yhat_sim)[1:3] <- c("a", "m", "p")

yhat_sim$m_repl <- yrepl_sim

ggplot(yhat_sim, aes(time, m, group=id)) + geom_line() + geom_point(aes(y=m_repl)) + facet_wrap(~id)


## now fit the simulated data set in parallel and serially, then compare


fake_data_serial_rk45 <- modifyList(data_serial, list(observed=1, yobs=yrepl_sim, solver=0))
fake_data_serial_bdf <- modifyList(data_serial, list(observed=1, yobs=yrepl_sim, solver=1))
fake_data_par_rk45 <- modifyList(fake_data_serial_rk45, list(parallel=1, solver=0))
fake_data_par_bdf <- modifyList(fake_data_par_rk45, list(parallel=1, solver=1))
fake_data_par_am <- modifyList(fake_data_par_rk45, list(parallel=1, solver=2))

fake_data_par_am$pbar

control <- list(adapt_delta=0.8, stepsize=1)

##Sys.setenv(OMP_NUM_THREADS = "4")
Sys.getenv("OMP_NUM_THREADS")


system.time(fake_run_par_rk45 <- sampling(sm, data=fake_data_par_rk45, warmup=250, iter=400, chains=1, init=list(truth), seed=12, refresh=1,
                                     control=control))

system.time(fake_run_par_bdf <- sampling(sm, data=fake_data_par_bdf, warmup=250, iter=400, chains=1, init=list(truth), seed=12, refresh=1,
                                     control=control))

system.time(fake_run_par_am <- sampling(sm, data=fake_data_par_am, warmup=250, iter=400, chains=1, init=list(truth), seed=12, refresh=1,
                                     control=control))

## compare with serial code
system.time(fake_run_serial_rk45 <- sampling(sm, data=fake_data_serial_rk45, warmup=250, iter=400, chains=1, init=list(truth), seed=12, refresh=1,
                                        control=control))

system.time(fake_run_serial_bdf <- sampling(sm, data=fake_data_serial_bdf, warmup=250, iter=400, chains=1, init=list(truth), seed=12, refresh=1,
                                        control=control))


## now look at sampler parameters and sample quality
colMeans(get_sampler_params(fake_run_par_rk45, inc_warmup=FALSE)[[1]])

print(fake_run_par_rk45, c("theta_v", "dose0_v"))

colMeans(get_sampler_params(fake_run_par_bdf, inc_warmup=FALSE)[[1]])

print(fake_run_par_bdf, c("theta_v", "dose0_v"))

colMeans(get_sampler_params(fake_run_par_am, inc_warmup=FALSE)[[1]])

print(fake_run_par_am, c("theta_v", "dose0_v"))

colMeans(get_sampler_params(fake_run_serial_rk45, inc_warmup=FALSE)[[1]])

print(fake_run_serial_rk45, c("theta_v", "dose0_v"))

colMeans(get_sampler_params(fake_run_serial_bdf, inc_warmup=FALSE)[[1]])

print(fake_run_serial_bdf, c("theta_v", "dose0_v"))

pars <- c("theta_v", "dose0_v", "sigma_y")

post_par <- as.matrix(fake_run_par_bdf, pars=pars)
post_serial <- as.matrix(fake_run_serial_rk45, pars=pars)

all(post_par == post_serial)
sum(post_par -  post_serial)

## compare estimates vs the truth
library(bayesplot)

true <- as.vector(unlist(truth))

bias_par <- sweep(post_par, 2, true)

mcmc_recover_intervals(bias_par, rep(0, length(true)))

bias_serial <- sweep(post_serial, 2, true)

mcmc_recover_intervals(bias_serial, rep(0, length(true)))
