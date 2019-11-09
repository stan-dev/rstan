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

## ----more-knitr-ops, include=FALSE---------------------------------------
knitr::opts_chunk$set(
  cache=TRUE,
  message=FALSE, 
  warning=FALSE
)

## ----lpdf, eval=FALSE----------------------------------------------------
#  /**
#   * Normal log-pdf for spatially lagged responses
#   *
#   * @param y Vector of response values.
#   * @param mu Mean parameter vector.
#   * @param sigma Positive scalar residual standard deviation.
#   * @param rho Positive scalar autoregressive parameter.
#   * @param W Spatial weight matrix.
#   *
#   * @return A scalar to be added to the log posterior.
#   */
#  real normal_lagsar_lpdf(vector y, vector mu, real sigma,
#                          real rho, matrix W) {
#    int N = rows(y);
#    real inv_sigma2 = 1 / square(sigma);
#    matrix[N, N] W_tilde = -rho * W;
#    vector[N] half_pred;
#  
#    for (n in 1:N) W_tilde[n,n] += 1;
#  
#    half_pred = W_tilde * (y - mdivide_left(W_tilde, mu));
#  
#    return 0.5 * log_determinant(crossprod(W_tilde) * inv_sigma2) -
#           0.5 * dot_self(half_pred) * inv_sigma2;
#  }

## ----setup, cache=FALSE--------------------------------------------------
library("loo")
library("brms")
library("bayesplot")
library("ggplot2")
color_scheme_set("brightblue")
theme_set(theme_default())


SEED <- 10001 
set.seed(SEED) # only sets seed for R (seed for Stan set later)

# loads COL.OLD data frame and COL.nb neighbor list
data(oldcol, package = "spdep") 

## ----data----------------------------------------------------------------
str(COL.OLD[, c("CRIME", "HOVAL", "INC")])

## ----fit, results="hide"-------------------------------------------------
fit <- brm(
  CRIME ~ INC + HOVAL, 
  data = COL.OLD,
  autocor = cor_lagsar(COL.nb),
  chains = 4,
  seed = SEED
)

## ----plot-lagsar, message=FALSE------------------------------------------
lagsar <- as.matrix(fit, pars = "lagsar")
estimates <- quantile(lagsar, probs = c(0.25, 0.5, 0.75))
mcmc_hist(lagsar) + 
  vline_at(estimates, linetype = 2, size = 1) +
  ggtitle("lagsar: posterior median and 50% central interval")

## ----approx--------------------------------------------------------------
posterior <- as.data.frame(fit)
y <- fit$data$CRIME
N <- length(y)
S <- nrow(posterior)
loglik <- yloo <- sdloo <- matrix(nrow = S, ncol = N)

for (s in 1:S) {
  p <- posterior[s, ]
  eta <- p$b_Intercept + p$b_INC * fit$data$INC + p$b_HOVAL * fit$data$HOVAL
  W_tilde <- diag(N) - p$lagsar * fit$autocor$W
  Cinv <- t(W_tilde) %*% W_tilde / p$sigma^2
  g <- Cinv %*% (y - solve(W_tilde, eta))
  cbar <- diag(Cinv)
  yloo[s, ] <- y - g / cbar
  sdloo[s, ] <- sqrt(1 / cbar)
  loglik[s, ] <- dnorm(y, yloo[s, ], sdloo[s, ], log = TRUE)
}

# use loo for psis smoothing
log_ratios <- -loglik
psis_result <- psis(log_ratios)

## ----plot, cache = FALSE-------------------------------------------------
plot(psis_result, label_points = TRUE)

## ----checklast, cache = FALSE--------------------------------------------
yloo_sub <- yloo[S, ]
sdloo_sub <- sdloo[S, ]
df <- data.frame(
  y = y, 
  yloo = yloo_sub,
  ymin = yloo_sub - sdloo_sub * 2,
  ymax = yloo_sub + sdloo_sub * 2
)
ggplot(data=df, aes(x = y, y = yloo, ymin = ymin, ymax = ymax)) +
  geom_errorbar(
    width = 1, 
    color = "skyblue3", 
    position = position_jitter(width = 0.25)
  ) +
  geom_abline(color = "gray30", size = 1.2) +
  geom_point()

## ----psisloo-------------------------------------------------------------
(psis_loo <- loo(loglik))

## ----fit_dummy, cache = TRUE---------------------------------------------
# see help("mi", "brms") for details on the mi() usage
fit_dummy <- brm(
  CRIME | mi() ~ INC + HOVAL, 
  data = COL.OLD,
  autocor = cor_lagsar(COL.nb), 
  chains = 0
)

## ----exact-loo-cv, results="hide", message=FALSE, warning=FALSE, cache = TRUE----
S <- 500
res <- vector("list", N)
loglik <- matrix(nrow = S, ncol = N)
for (i in seq_len(N)) {
  dat_mi <- COL.OLD
  dat_mi$CRIME[i] <- NA
  fit_i <- update(fit_dummy, newdata = dat_mi, 
                  # just for vignette
                  chains = 1, iter = S * 2)
  posterior <- as.data.frame(fit_i)
  yloo <- sdloo <- rep(NA, S)
  for (s in seq_len(S)) {
    p <- posterior[s, ]
    y_miss_i <- y
    y_miss_i[i] <- p$Ymi
    eta <- p$b_Intercept + p$b_INC * fit_i$data$INC + p$b_HOVAL * fit_i$data$HOVAL
    W_tilde <- diag(N) - p$lagsar * fit_i$autocor$W
    Cinv <- t(W_tilde) %*% W_tilde / p$sigma^2
    g <- Cinv %*% (y_miss_i - solve(W_tilde, eta))
    cbar <- diag(Cinv);
    yloo[s] <- y_miss_i[i] - g[i] / cbar[i]
    sdloo[s] <- sqrt(1 / cbar[i])
    loglik[s, i] <- dnorm(y[i], yloo[s], sdloo[s], log = TRUE)
  }
  ypred <- rnorm(S, yloo, sdloo)
  res[[i]] <- data.frame(y = c(posterior$Ymi, ypred))
  res[[i]]$type <- rep(c("pp", "loo"), each = S)
  res[[i]]$obs <- i
}
res <- do.call(rbind, res)

## ----yplots, cache = FALSE, fig.width=10, out.width="95%", fig.asp = 0.3----
res_sub <- res[res$obs %in% 1:4, ]
ggplot(res_sub, aes(y, fill = type)) +
  geom_density(alpha = 0.6) +
  facet_wrap("obs", scales = "fixed", ncol = 4)

## ----loo_exact, cache=FALSE----------------------------------------------
log_mean_exp <- function(x) {
  # more stable than log(mean(exp(x)))
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x))) - log(length(x))
}
exact_elpds <- apply(loglik, 2, log_mean_exp)
exact_elpd <- sum(exact_elpds)
round(exact_elpd, 1)

## ----compare, fig.height=5-----------------------------------------------
df <- data.frame(
  approx_elpd = psis_loo$pointwise[, "elpd_loo"],
  exact_elpd = exact_elpds
)
ggplot(df, aes(x = approx_elpd, y = exact_elpd)) +
  geom_abline(color = "gray30") +
  geom_point(size = 2) +
  geom_point(data = df[4, ], size = 3, color = "red3") +
  xlab("Approximate elpds") +
  ylab("Exact elpds") +
  coord_fixed(xlim = c(-16, -3), ylim = c(-16, -3))

## ----pt4-----------------------------------------------------------------
without_pt_4 <- c(
  approx = sum(psis_loo$pointwise[-4, "elpd_loo"]),
  exact = sum(exact_elpds[-4])  
)
round(without_pt_4, 1)

## ----brms-stan-code, eval=FALSE------------------------------------------
#  // generated with brms 2.2.0
#  functions {
#  /**
#   * Normal log-pdf for spatially lagged responses
#   *
#   * @param y Vector of response values.
#   * @param mu Mean parameter vector.
#   * @param sigma Positive scalar residual standard deviation.
#   * @param rho Positive scalar autoregressive parameter.
#   * @param W Spatial weight matrix.
#   *
#   * @return A scalar to be added to the log posterior.
#   */
#    real normal_lagsar_lpdf(vector y, vector mu, real sigma,
#                            real rho, matrix W) {
#      int N = rows(y);
#      real inv_sigma2 = 1 / square(sigma);
#      matrix[N, N] W_tilde = -rho * W;
#      vector[N] half_pred;
#      for (n in 1:N) W_tilde[n, n] += 1;
#      half_pred = W_tilde * (y - mdivide_left(W_tilde, mu));
#      return 0.5 * log_determinant(crossprod(W_tilde) * inv_sigma2) -
#             0.5 * dot_self(half_pred) * inv_sigma2;
#    }
#  }
#  data {
#    int<lower=1> N;  // total number of observations
#    vector[N] Y;  // response variable
#    int<lower=0> Nmi;  // number of missings
#    int<lower=1> Jmi[Nmi];  // positions of missings
#    int<lower=1> K;  // number of population-level effects
#    matrix[N, K] X;  // population-level design matrix
#    matrix[N, N] W;  // spatial weight matrix
#    int prior_only;  // should the likelihood be ignored?
#  }
#  transformed data {
#    int Kc = K - 1;
#    matrix[N, K - 1] Xc;  // centered version of X
#    vector[K - 1] means_X;  // column means of X before centering
#    for (i in 2:K) {
#      means_X[i - 1] = mean(X[, i]);
#      Xc[, i - 1] = X[, i] - means_X[i - 1];
#    }
#  }
#  parameters {
#    vector[Nmi] Ymi;  // estimated missings
#    vector[Kc] b;  // population-level effects
#    real temp_Intercept;  // temporary intercept
#    real<lower=0> sigma;  // residual SD
#    real<lower=0,upper=1> lagsar;  // SAR parameter
#  }
#  transformed parameters {
#  }
#  model {
#    vector[N] Yl = Y;
#    vector[N] mu = Xc * b + temp_Intercept;
#    Yl[Jmi] = Ymi;
#    // priors including all constants
#    target += student_t_lpdf(temp_Intercept | 3, 34, 17);
#    target += student_t_lpdf(sigma | 3, 0, 17)
#      - 1 * student_t_lccdf(0 | 3, 0, 17);
#    // likelihood including all constants
#    if (!prior_only) {
#      target += normal_lagsar_lpdf(Yl | mu, sigma, lagsar, W);
#    }
#  }
#  generated quantities {
#    // actual population-level intercept
#    real b_Intercept = temp_Intercept - dot_product(means_X, b);
#  }

