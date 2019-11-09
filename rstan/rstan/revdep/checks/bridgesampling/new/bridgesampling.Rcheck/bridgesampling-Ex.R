pkgname <- "bridgesampling"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bridgesampling')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bridge_sampler")
### * bridge_sampler

flush(stderr()); flush(stdout())

### Name: bridge_sampler
### Title: Log Marginal Likelihood via Bridge Sampling
### Aliases: bridge_sampler bridge_sampler.stanfit bridge_sampler.mcmc.list
###   bridge_sampler.mcmc bridge_sampler.matrix bridge_sampler.stanreg
###   bridge_sampler.rjags bridge_sampler.runjags
###   bridge_sampler.MCMC_refClass

### ** Examples

## ------------------------------------------------------------------------
## Example 1: Estimating the Normalizing Constant of a Two-Dimensional
##            Standard Normal Distribution
## ------------------------------------------------------------------------

library(bridgesampling)
library(mvtnorm)

samples <- rmvnorm(1e4, mean = rep(0, 2), sigma = diag(2))
colnames(samples) <- c("x1", "x2")
log_density <- function(samples.row, data) {
  -.5*t(samples.row) %*% samples.row
}

lb <- rep(-Inf, 2)
ub <- rep(Inf, 2)
names(lb) <- names(ub) <- colnames(samples)
bridge_result <- bridge_sampler(samples = samples, log_posterior = log_density,
                                data = NULL, lb = lb, ub = ub, silent = TRUE)

# compare to analytical value
analytical <- log(2*pi)
print(cbind(bridge_result$logml, analytical))

## Not run: 
##D 
##D ## ------------------------------------------------------------------------
##D ## Example 2: Hierarchical Normal Model
##D ## ------------------------------------------------------------------------
##D 
##D # for a full description of the example, see
##D vignette("bridgesampling_example_jags")
##D 
##D library(R2jags)
##D 
##D ### generate data ###
##D 
##D set.seed(12345)
##D 
##D mu <- 0
##D tau2 <- 0.5
##D sigma2 <- 1
##D 
##D n <- 20
##D theta <- rnorm(n, mu, sqrt(tau2))
##D y <- rnorm(n, theta, sqrt(sigma2))
##D 
##D 
##D ### set prior parameters
##D alpha <- 1
##D beta <- 1
##D mu0 <- 0
##D tau20 <- 1
##D 
##D ### functions to get posterior samples ###
##D 
##D ### H0: mu = 0
##D 
##D getSamplesModelH0 <- function(data, niter = 52000, nburnin = 2000, nchains = 3) {
##D 
##D   model <- "
##D     model {
##D       for (i in 1:n) {
##D         theta[i] ~ dnorm(0, invTau2)
##D           y[i] ~ dnorm(theta[i], 1/sigma2)
##D       }
##D       invTau2 ~ dgamma(alpha, beta)
##D       tau2 <- 1/invTau2
##D     }"
##D 
##D   s <- jags(data, parameters.to.save = c("theta", "invTau2"),
##D             model.file = textConnection(model),
##D             n.chains = nchains, n.iter = niter,
##D             n.burnin = nburnin, n.thin = 1)
##D 
##D   return(s)
##D 
##D }
##D 
##D ### H1: mu != 0
##D 
##D getSamplesModelH1 <- function(data, niter = 52000, nburnin = 2000,
##D                               nchains = 3) {
##D 
##D   model <- "
##D     model {
##D       for (i in 1:n) {
##D         theta[i] ~ dnorm(mu, invTau2)
##D         y[i] ~ dnorm(theta[i], 1/sigma2)
##D       }
##D       mu ~ dnorm(mu0, 1/tau20)
##D       invTau2 ~ dgamma(alpha, beta)
##D       tau2 <- 1/invTau2
##D     }"
##D 
##D   s <- jags(data, parameters.to.save = c("theta", "mu", "invTau2"),
##D             model.file = textConnection(model),
##D             n.chains = nchains, n.iter = niter,
##D             n.burnin = nburnin, n.thin = 1)
##D 
##D   return(s)
##D 
##D }
##D 
##D ### get posterior samples ###
##D 
##D # create data lists for Jags
##D data_H0 <- list(y = y, n = length(y), alpha = alpha, beta = beta, sigma2 = sigma2)
##D data_H1 <- list(y = y, n = length(y), mu0 = mu0, tau20 = tau20, alpha = alpha,
##D                 beta = beta, sigma2 = sigma2)
##D 
##D # fit models
##D samples_H0 <- getSamplesModelH0(data_H0)
##D samples_H1 <- getSamplesModelH1(data_H1)
##D 
##D 
##D ### functions for evaluating the unnormalized posteriors on log scale ###
##D log_posterior_H0 <- function(samples.row, data) {
##D 
##D   mu <- 0
##D   invTau2 <- samples.row[[ "invTau2" ]]
##D   theta <- samples.row[ paste0("theta[", seq_along(data$y), "]") ]
##D 
##D   sum(dnorm(data$y, theta, data$sigma2, log = TRUE)) +
##D     sum(dnorm(theta, mu, 1/sqrt(invTau2), log = TRUE)) +
##D     dgamma(invTau2, data$alpha, data$beta, log = TRUE)
##D 
##D }
##D 
##D log_posterior_H1 <- function(samples.row, data) {
##D 
##D   mu <- samples.row[[ "mu" ]]
##D   invTau2 <- samples.row[[ "invTau2" ]]
##D   theta <- samples.row[ paste0("theta[", seq_along(data$y), "]") ]
##D 
##D   sum(dnorm(data$y, theta, data$sigma2, log = TRUE)) +
##D     sum(dnorm(theta, mu, 1/sqrt(invTau2), log = TRUE)) +
##D     dnorm(mu, data$mu0, sqrt(data$tau20), log = TRUE) +
##D     dgamma(invTau2, data$alpha, data$beta, log = TRUE)
##D 
##D }
##D 
##D # specify parameter bounds H0
##D cn <- colnames(samples_H0$BUGSoutput$sims.matrix)
##D cn <- cn[cn != "deviance"]
##D lb_H0 <- rep(-Inf, length(cn))
##D ub_H0 <- rep(Inf, length(cn))
##D names(lb_H0) <- names(ub_H0) <- cn
##D lb_H0[[ "invTau2" ]] <- 0
##D 
##D # specify parameter bounds H1
##D cn <- colnames(samples_H1$BUGSoutput$sims.matrix)
##D cn <- cn[cn != "deviance"]
##D lb_H1 <- rep(-Inf, length(cn))
##D ub_H1 <- rep(Inf, length(cn))
##D names(lb_H1) <- names(ub_H1) <- cn
##D lb_H1[[ "invTau2" ]] <- 0
##D 
##D 
##D # compute log marginal likelihood via bridge sampling for H0
##D H0.bridge <- bridge_sampler(samples = samples_H0, data = data_H0,
##D                             log_posterior = log_posterior_H0, lb = lb_H0,
##D                             ub = ub_H0, silent = TRUE)
##D print(H0.bridge)
##D 
##D # compute log marginal likelihood via bridge sampling for H1
##D H1.bridge <- bridge_sampler(samples = samples_H1, data = data_H1,
##D                             log_posterior = log_posterior_H1, lb = lb_H1,
##D                             ub = ub_H1, silent = TRUE)
##D print(H1.bridge)
##D 
##D # compute percentage error
##D print(error_measures(H0.bridge)$percentage)
##D print(error_measures(H1.bridge)$percentage)
##D 
##D # compute Bayes factor
##D BF01 <- bf(H0.bridge, H1.bridge)
##D print(BF01)
##D 
##D # compute posterior model probabilities (assuming equal prior model probabilities)
##D post1 <- post_prob(H0.bridge, H1.bridge)
##D print(post1)
##D 
##D # compute posterior model probabilities (using user-specified prior model probabilities)
##D post2 <- post_prob(H0.bridge, H1.bridge, prior_prob = c(.6, .4))
##D print(post2)
##D 
## End(Not run)

## Not run: 
##D 
##D ## ------------------------------------------------------------------------
##D ## Example 3: rstanarm
##D ## ------------------------------------------------------------------------
##D library(rstanarm)
##D 
##D # N.B.: remember to specify the diagnostic_file
##D 
##D fit_1 <- stan_glm(mpg ~ wt + qsec + am, data = mtcars,
##D                   chains = 2, cores = 2, iter = 5000,
##D                   diagnostic_file = file.path(tempdir(), "df.csv"))
##D bridge_1 <- bridge_sampler(fit_1)
##D fit_2 <- update(fit_1, formula = . ~ . + cyl)
##D bridge_2 <- bridge_sampler(fit_2, method = "warp3")
##D bf(bridge_1, bridge_2)
##D 
## End(Not run)




cleanEx()
nameEx("ier")
### * ier

flush(stderr()); flush(stdout())

### Name: ier
### Title: Standardized International Exchange Rate Changes from 1975 to
###   1986
### Aliases: ier
### Keywords: dataset

### ** Examples


## Not run: 
##D 
##D ################################################################################
##D # BAYESIAN FACTOR ANALYSIS (AS PROPOSED BY LOPES & WEST, 2004)
##D ################################################################################
##D 
##D library(bridgesampling)
##D library(rstan)
##D 
##D cores <- 4
##D options(mc.cores = cores)
##D 
##D data("ier")
##D 
##D #-------------------------------------------------------------------------------
##D # plot data
##D #-------------------------------------------------------------------------------
##D 
##D currency <- colnames(ier)
##D label <- c("US Dollar", "Canadian Dollar", "Yen", "Franc", "Lira", "Mark")
##D op <- par(mfrow = c(3, 2), mar = c(6, 6, 3, 3))
##D 
##D for (i in seq_along(currency)) {
##D   plot(ier[,currency[i]], type = "l", col = "darkblue",  axes = FALSE,
##D        ylim = c(-4, 4), ylab = "", xlab = "", lwd = 2)
##D   axis(1, at = 0:12*12, labels = 1975:1987, cex.axis = 1.7)
##D   axis(2, at = pretty(c(-4, 4)), las = 1, cex.axis = 1.7)
##D   mtext("Year", 1, cex = 1.5, line = 3.2)
##D   mtext("Exchange Rate Changes", 2, cex = 1.4, line = 3.2)
##D   mtext(label[i], 3, cex = 1.6, line = .1)
##D }
##D 
##D par(op)
##D 
##D #-------------------------------------------------------------------------------
##D # stan model
##D #-------------------------------------------------------------------------------
##D 
##D model_code <-
##D "data {
##D   int<lower=1> T; // number of observations
##D   int<lower=1> m; // number of variables
##D   int<lower=1> k; // number of factors
##D   matrix[T,m] Y;  // data matrix
##D }
##D transformed data {
##D   int<lower = 1> r;
##D   vector[m] zeros;
##D   r = m * k - k * (k - 1) / 2; // number of non-zero factor loadings
##D   zeros = rep_vector(0.0, m);
##D }
##D parameters {
##D   real beta_lower[r - k];  // lower-diagonal elements of beta
##D   real<lower = 0> beta_diag [k]; // diagonal elements of beta
##D   vector<lower = 0>[m] sigma2; // residual variances
##D }
##D transformed parameters {
##D   matrix[m,k] beta;
##D   cov_matrix[m] Omega;
##D   // construct lower-triangular factor loadings matrix
##D   {
##D     int index_lower = 1;
##D     for (j in 1:k) {
##D       for (i in 1:m) {
##D         if (i == j) {
##D           beta[j,j] = beta_diag[j];
##D         } else if (i >= j) {
##D           beta[i,j] = beta_lower[index_lower];
##D           index_lower = index_lower + 1;
##D         } else {
##D           beta[i,j] = 0.0;
##D         }
##D       }
##D     }
##D   }
##D   Omega = beta * beta' + diag_matrix(sigma2);
##D }
##D model {
##D   // priors
##D   target += normal_lpdf(beta_diag | 0, 1) - k * normal_lccdf(0 | 0, 1);
##D   target += normal_lpdf(beta_lower | 0, 1);
##D   target += inv_gamma_lpdf(sigma2 | 2.2 / 2.0, 0.1 / 2.0);
##D 
##D   // likelihood
##D   for(t in 1:T) {
##D     target += multi_normal_lpdf(Y[t] | zeros, Omega);
##D   }
##D }"
##D 
##D # compile model
##D model <- stan_model(model_code = model_code)
##D 
##D 
##D #-------------------------------------------------------------------------------
##D # fit models and compute log marginal likelihoods
##D #-------------------------------------------------------------------------------
##D 
##D # function for generating starting values
##D init_fun <- function(nchains, k, m) {
##D   r <- m * k - k * (k - 1) / 2
##D   out <- vector("list", nchains)
##D   for (i in seq_len(nchains)) {
##D     beta_lower <- array(runif(r - k, 0.05, 1), dim = r - k)
##D     beta_diag <- array(runif(k, .05, 1), dim = k)
##D     sigma2 <- array(runif(m, .05, 1.5), dim = m)
##D     out[[i]] <- list(beta_lower = beta_lower,
##D                      beta_diag = beta_diag,
##D                      sigma2 = sigma2)
##D   }
##D   return(out)
##D }
##D 
##D set.seed(1)
##D stanfit <- bridge <- vector("list", 3)
##D for (k in 1:3) {
##D   stanfit[[k]] <- sampling(model,
##D                            data = list(Y = ier, T = nrow(ier),
##D                                        m = ncol(ier), k = k),
##D                            iter = 11000, warmup = 1000, chains = 4,
##D                            init = init_fun(nchains = 4, k = k, m = ncol(ier)),
##D                            cores = cores, seed = 1)
##D   bridge[[k]] <- bridge_sampler(stanfit[[k]], method = "warp3",
##D                                 repetitions = 10, cores = cores)
##D }
##D 
##D # example output
##D summary(bridge[[2]])
##D 
##D #-------------------------------------------------------------------------------
##D # compute posterior model probabilities
##D #-------------------------------------------------------------------------------
##D 
##D pp <- post_prob(bridge[[1]], bridge[[2]], bridge[[3]],
##D           model_names = c("k = 1", "k = 2", "k = 3"))
##D pp
##D 
##D op <- par(mar = c(6, 6, 3, 3))
##D boxplot(pp, axes = FALSE,
##D      ylim = c(0, 1), ylab = "",
##D      xlab = "")
##D axis(1, at = 1:3, labels = colnames(pp), cex.axis = 1.7)
##D axis(2, cex.axis = 1.1)
##D mtext("Posterior Model Probability", 2, cex = 1.5, line = 3.2)
##D mtext("Number of Factors", 1, cex = 1.4, line = 3.2)
##D par(op)
##D 
## End(Not run)



cleanEx()
nameEx("post_prob")
### * post_prob

flush(stderr()); flush(stdout())

### Name: post_prob
### Title: Posterior Model Probabilities from Marginal Likelihoods
### Aliases: post_prob post_prob.bridge post_prob.bridge_list
###   post_prob.default

### ** Examples


H0 <- structure(list(logml = -20.8084543022433, niter = 4, method = "normal"),
                .Names = c("logml", "niter", "method"), class = "bridge")
H1 <- structure(list(logml = -17.9623077558729, niter = 4, method = "normal"),
                .Names = c("logml", "niter", "method"), class = "bridge")
H2 <- structure(list(logml = -19, niter = 4, method = "normal"),
                .Names = c("logml", "niter", "method"), class = "bridge")


post_prob(H0, H1, H2)
post_prob(H1, H0)

## all produce the same (only names differ):
post_prob(H0, H1, H2)
post_prob(H0$logml, H1$logml, H2$logml)
post_prob(c(H0$logml, H1$logml, H2$logml))
post_prob(H0$logml, c(H1$logml, H2$logml))
post_prob(H0$logml, c(H1$logml, H2$logml), model_names = c("H0", "H1", "H2"))


### with bridge list elements:
H0L <- structure(list(logml = c(-20.8088381186739, -20.8072772698116,
-20.808454454621, -20.8083419072281, -20.8087870541247, -20.8084887398113,
-20.8086023582344, -20.8079083169745, -20.8083048489095, -20.8090050811436
), niter = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4), method = "normal",
    repetitions = 10), .Names = c("logml", "niter", "method",
"repetitions"), class = "bridge_list")

H1L <- structure(list(logml = c(-17.961665507006, -17.9611290723151,
-17.9607509604499, -17.9608629535992, -17.9602093576442, -17.9600223300432,
-17.9610157118017, -17.9615557696561, -17.9608437034849, -17.9606743200309
), niter = c(4, 4, 4, 4, 4, 4, 4, 4, 3, 4), method = "normal",
    repetitions = 10), .Names = c("logml", "niter", "method",
"repetitions"), class = "bridge_list")

post_prob(H1L, H0L)
post_prob(H1L, H0L, H0) # last element recycled with warning.




cleanEx()
nameEx("turtles")
### * turtles

flush(stderr()); flush(stdout())

### Name: turtles
### Title: Turtles Data from Janzen, Tucker, and Paukstis (2000)
### Aliases: turtles
### Keywords: dataset

### ** Examples


## Not run: 
##D 
##D ################################################################################
##D # BAYESIAN GENERALIZED LINEAR MIXED MODEL (PROBIT REGRESSION)
##D ################################################################################
##D 
##D library(bridgesampling)
##D library(rstan)
##D 
##D data("turtles")
##D 
##D #-------------------------------------------------------------------------------
##D # plot data
##D #-------------------------------------------------------------------------------
##D 
##D # reproduce Figure 1 from Sinharay & Stern (2005)
##D xticks <- pretty(turtles$clutch)
##D yticks <- pretty(turtles$x)
##D 
##D plot(1, type = "n", axes = FALSE, ylab = "", xlab = "", xlim = range(xticks),
##D      ylim =  range(yticks))
##D points(turtles$clutch, turtles$x, pch = ifelse(turtles$y, 21, 4), cex = 1.3,
##D        col = ifelse(turtles$y, "black", "darkred"), bg = "grey", lwd = 1.3)
##D axis(1, cex.axis = 1.4)
##D mtext("Clutch Identifier", side = 1, line = 2.9, cex = 1.8)
##D axis(2, las = 1, cex.axis = 1.4)
##D mtext("Birth Weight (Grams)", side = 2, line = 2.6, cex = 1.8)
##D 
##D #-------------------------------------------------------------------------------
##D # Analysis: Natural Selection Study (compute same BF as Sinharay & Stern, 2005)
##D #-------------------------------------------------------------------------------
##D 
##D ### H0 (model without random intercepts) ###
##D H0_code <-
##D "data {
##D   int<lower = 1> N;
##D   int<lower = 0, upper = 1> y[N];
##D   real<lower = 0> x[N];
##D }
##D parameters {
##D   real alpha0_raw;
##D   real alpha1_raw;
##D }
##D transformed parameters {
##D   real alpha0 = sqrt(10.0) * alpha0_raw;
##D   real alpha1 = sqrt(10.0) * alpha1_raw;
##D }
##D model {
##D   // priors
##D   target += normal_lpdf(alpha0_raw | 0, 1);
##D   target += normal_lpdf(alpha1_raw | 0, 1);
##D 
##D   // likelihood
##D   for (i in 1:N) {
##D     target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1 * x[i]));
##D   }
##D }"
##D 
##D ### H1 (model with random intercepts) ###
##D H1_code <-
##D "data {
##D   int<lower = 1> N;
##D   int<lower = 0, upper = 1> y[N];
##D   real<lower = 0> x[N];
##D   int<lower = 1> C;
##D   int<lower = 1, upper = C> clutch[N];
##D }
##D parameters {
##D   real alpha0_raw;
##D   real alpha1_raw;
##D   vector[C] b_raw;
##D   real<lower = 0> sigma2;
##D }
##D transformed parameters {
##D   vector[C] b;
##D   real<lower = 0> sigma = sqrt(sigma2);
##D   real alpha0 = sqrt(10.0) * alpha0_raw;
##D   real alpha1 = sqrt(10.0) * alpha1_raw;
##D   b = sigma * b_raw;
##D }
##D model {
##D   // priors
##D   target += - 2 * log(1 + sigma2); // p(sigma2) = 1 / (1 + sigma2) ^ 2
##D   target += normal_lpdf(alpha0_raw | 0, 1);
##D   target += normal_lpdf(alpha1_raw | 0, 1);
##D 
##D   // random effects
##D   target += normal_lpdf(b_raw | 0, 1);
##D 
##D   // likelihood
##D   for (i in 1:N) {
##D     target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1 * x[i] + b[clutch[i]]));
##D   }
##D }"
##D 
##D set.seed(1)
##D ### fit models ###
##D stanfit_H0 <- stan(model_code = H0_code,
##D                    data = list(y = turtles$y, x = turtles$x, N = nrow(turtles)),
##D                    iter = 15500, warmup = 500, chains = 4, seed = 1)
##D stanfit_H1 <- stan(model_code = H1_code,
##D                    data = list(y = turtles$y, x = turtles$x, N = nrow(turtles),
##D                                C = max(turtles$clutch), clutch = turtles$clutch),
##D                    iter = 15500, warmup = 500, chains = 4, seed = 1)
##D 
##D set.seed(1)
##D ### compute (log) marginal likelihoods ###
##D bridge_H0 <- bridge_sampler(stanfit_H0)
##D bridge_H1 <- bridge_sampler(stanfit_H1)
##D 
##D ### compute approximate percentage errors ###
##D error_measures(bridge_H0)$percentage
##D error_measures(bridge_H1)$percentage
##D 
##D ### summary ###
##D summary(bridge_H0)
##D summary(bridge_H1)
##D 
##D ### compute Bayes factor ("true" value: BF01 = 1.273) ###
##D bf(bridge_H0, bridge_H1)
##D 
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
