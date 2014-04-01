.setUp <- function() { }
test_output_csv_and_extract <- function() {
  csv_fname <- 'teststanfit.csv'
  model_code <- "
    transformed data {
      int n;
      n <- 0;
    }
    parameters { 
      real y[2];
      real zeroleny[n];
    } 
    transformed parameters {
      real y2[2, 2];
      y2[1, 1] <- y[1]; 
      y2[1, 2] <- -y[1]; 
      y2[2, 1] <- -y[2]; 
      y2[2, 2] <- y[2]; 
    } 
    model {
      y ~ normal(0, 1);
    } 
    generated quantities {
      real y3[3, 2, 3];
      y3[1,1,1] <- 1;   y3[1,1,2] <- 2;   y3[1,1,3] <- 3;  
      y3[1,2,1] <- 4;   y3[1,2,2] <- 5;   y3[1,2,3] <- 6;
      y3[2,1,1] <- 7;   y3[2,1,2] <- 8;   y3[2,1,3] <- 9;  
      y3[2,2,1] <- 10;  y3[2,2,2] <- 11;  y3[2,2,3] <- 12;
      y3[3,1,1] <- 13;  y3[3,1,2] <- 14;  y3[3,1,3] <- 15;  
      y3[3,2,1] <- 16;  y3[3,2,2] <- 17;  y3[3,2,3] <- 18;
    } 
  "

  ## Disable the zero length vector as there is a bug with relist in R (< 3.0.3).
  ## See: https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=15499
  ## https://github.com/stan-dev/rstan/issues/51
  rversion_date <- 
    as.Date(paste(R.version$year, R.version$month, R.version$day, sep = '-')) 
  rversion_dat_302 <- as.Date("2014-03-06")
  if (rversion_date < rversion_dat_302) {
    model_code <- gsub('.*zeroleny.*', '', model_code, perl = TRUE)
  } 
  ## 
 
  fit <- stan(model_code = model_code, 
              iter = 100, chains = 1, thin = 3, 
              sample_file = csv_fname)

  pmean <- get_posterior_mean(fit)[,1]
  pmean2 <- summary(fit)$c_summary[,'mean',1]
  checkEquals(pmean, pmean2, tolerance = 0.00001)

  checkTrue(file.exists(csv_fname))
  d <- read.csv(file = csv_fname, comment.char = '#',
                header = TRUE)
  y_names <- paste0("y.", 1:2)
  cat2 <- function(a, b) paste0(a, ".", b)
  y2_names <- paste0("y2.", t(outer(1:2, 1:2, cat2)))
  y3_names <- paste0("y3.", t(outer(as.vector(t(outer(1:3, 1:2, cat2))), 1:3, cat2)))

  dfit <- as.data.frame(fit)
  iter1 <- unlist(d[1,])
  names(iter1) <-  rstan:::sqrfnames_to_dotfnames(names(iter1))
  checkEquals(iter1["y.1"], iter1["y2.1.1"], checkNames = FALSE)
  checkEquals(iter1["y.1"], -iter1["y2.1.2"], checkNames = FALSE)
  checkEquals(iter1["y.2"], iter1["y2.2.2"], checkNames = FALSE)
  checkEquals(iter1["y.2"], -iter1["y2.2.1"], checkNames = FALSE)
  checkEquals(iter1[y3_names], 1:18, checkNames = FALSE)

  # FIXME, uncomment the follwing line
  # checkEquals(colnames(d), c("lp__", "treedepth__", "stepsize__", y_names, y2_names, y3_names))
  iter1 <- unlist(d[1,])
  # to check if the order is also row-major
  checkEquals(iter1["y.1"], iter1["y2.1.1"], checkNames = FALSE)
  checkEquals(iter1["y.1"], -iter1["y2.1.2"], checkNames = FALSE)
  checkEquals(iter1["y.2"], iter1["y2.2.2"], checkNames = FALSE)
  checkEquals(iter1["y.2"], -iter1["y2.2.1"], checkNames = FALSE)
  checkEquals(iter1[y3_names], 1:18, checkNames = FALSE)

  fit2 <- read_stan_csv(csv_fname)
  checkEquals(rstan:::is_sf_valid(fit), TRUE)
  checkEquals(rstan:::is_sf_valid(fit2), FALSE)
  fit2_m <- as.vector(get_posterior_mean(fit2))
  fit2_m2 <- summary(fit2)$summary[,"mean"]
  checkEquals(fit2_m, fit2_m2, checkNames = FALSE)
  unlink(csv_fname)

  # test extract 
  fitb <- stan(fit = fit, iter = 105, warmup = 100, chains = 3, thin = 3)
  e1 <- extract(fitb)
  checkEquals(e1$y3[1,1,1,1], 1, checkNames = FALSE)
  checkEquals(e1$y3[1,2,2,2], 11, checkNames = FALSE)
  checkEquals(e1$y3[1,2,1,2], 8, checkNames = FALSE)
  checkEquals(e1$y3[1,3,1,2], 14, checkNames = FALSE)
} 

test_domain_error_exception <- function() { 
  code <- '
  parameters {
    real x;
  }
  model {
    x ~ normal(0, -1);
  }
  '
  sm <- stan_model(model_code = code)
  mod <- sm@dso@.CXXDSOMISC$module
  model_cppname <- sm@model_cpp$model_cppname
  sf_mod <- eval(call("$", mod, paste0('stan_fit4', model_cppname)))
  dat <- list()
  sampler <- new(sf_mod, dat, sm@dso@.CXXDSOMISC$cxxfun)
  args <- list(init = "0", iter = 1)
  checkException(s <- sampler$call_sampler(args))
  checkTrue(grepl('.*Domain error during initialization with 0.*', geterrmessage()))
} 

test_init_zero_exception_inf_lp <- function() {
  code <- '
  parameters {
    real x;
  }
  model {
    lp__ <- 1 / x;
  }
  '
  sm <- stan_model(model_code = code)
  mod <- sm@dso@.CXXDSOMISC$module
  model_cppname <- sm@model_cpp$model_cppname
  sf_mod <- eval(call("$", mod, paste0('stan_fit4', model_cppname)))
  dat <- list()
  sampler <- new(sf_mod, dat, sm@dso@.CXXDSOMISC$cxxfun)
  args <- list(init = "0", iter = 1)
  checkException(s <- sampler$call_sampler(args))
  checkTrue(grepl('.*vanishing density.*', geterrmessage()))
} 

test_init_zero_exception_inf_grad <- function() {
  code <- '
  parameters {
    real x;
  }
  model {
    lp__ <- 1 / log(x);
  }
  '
  sm <- stan_model(model_code = code)
  mod <- sm@dso@.CXXDSOMISC$module
  model_cppname <- sm@model_cpp$model_cppname
  sf_mod <- eval(call("$", mod, paste0('stan_fit4', model_cppname)))
  dat <- list()
  sampler <- new(sf_mod, dat, sm@dso@.CXXDSOMISC$cxxfun)
  args <- list(init = "0", iter = 1)
  checkException(s <- sampler$call_sampler(args))
  checkTrue(grepl('.*divergent gradient.*', geterrmessage()))
}

test_grad_log <- function() {
  y <- c(0.70,  -0.16,  0.77, -1.37, -1.99,  1.35, 0.08, 
         0.02,  -1.48, -0.08,  0.34,  0.03, -0.42, 0.87, 
         -1.36,  1.43,  0.80, -0.48, -1.61, -1.27)

  code <- '
  data {
    real y[20];
  } 
  parameters {
    real mu;
    real<lower=0> sigma;
  } 
  model {
    y ~ normal(mu, sigma);
  } 
  '
  log_prob_fun <- function(mu, log_sigma, adjust = TRUE) {
    sigma <- exp(log_sigma)
    lp <- -sum((y - mu)^2) / (2 * (sigma^2)) - length(y) * log(sigma) 
    if (adjust) lp <- lp + log(sigma)
    lp
  } 
  log_prob_grad_fun <- function(mu, log_sigma, adjust = TRUE) {
    sigma <- exp(log_sigma)
    g_lsigma <- sum((y - mu)^2) * sigma^(-2) - length(y) 
    if (adjust) g_lsigma <- g_lsigma + 1
    g_mu <- sum(y - mu) * sigma^(-2)
    c(g_mu, g_lsigma)
  } 
  sf <- stan(model_code = code, data = list(y = y), iter = 200)
  mu <- 0.1; sigma <- 2;
  checkEquals(log_prob(sf, unconstrain_pars(sf, list(mu = mu, sigma = sigma))), 
              log_prob_fun(mu, log(sigma)), checkNames = FALSE)
  checkEquals(log_prob(sf, unconstrain_pars(sf, list(mu = mu, sigma = sigma)), FALSE), 
              log_prob_fun(mu, log(sigma), adjust = FALSE), checkNames = FALSE)
  lp1 <- log_prob(sf, unconstrain_pars(sf, list(mu = mu, sigma = sigma)), FALSE, TRUE)
  checkEquals(attr(lp1, 'gradient'), log_prob_grad_fun(mu, log(sigma), FALSE))
  g1 <- grad_log_prob(sf, unconstrain_pars(sf, list(mu = mu, sigma = sigma)), FALSE)
  checkEquals(attr(g1, 'log_prob'), log_prob_fun(mu, log(sigma), adjust = FALSE))
  attributes(g1) <- NULL
  checkEquals(g1, log_prob_grad_fun(mu, log(sigma), adjust = FALSE))
}

test_specify_args <- function() {
  y <- c(0.70,  -0.16,  0.77, -1.37, -1.99,  1.35, 0.08, 
         0.02,  -1.48, -0.08,  0.34,  0.03, -0.42, 0.87, 
         -1.36,  1.43,  0.80, -0.48, -1.61, -1.27)

  code <- '
  data {
    real y[20];
  } 
  parameters {
    real mu;
    real<lower=0> sigma;
  } 
  model {
    y ~ normal(mu, sigma);
  } 
  '
  stepsize0 <- 0.15
  sf <- stan(model_code = code, data = list(y = y), iter = 200, 
             control = list(adapt_engaged = FALSE, stepsize = stepsize0))
  checkEquals(attr(sf@sim$samples[[1]],"sampler_params")$stepsize__[1], stepsize0)

  sf2 <- stan(fit = sf, iter = 20, algorithm = 'HMC', data = list(y = y),
             control = list(adapt_engaged = FALSE, stepsize = stepsize0))
  checkEquals(attr(sf2@sim$samples[[1]],"sampler_params")$stepsize__[1], stepsize0)

  sf3 <- stan(fit = sf, iter = 1, data = list(y = y), init = 0, chains = 1)
  i_u <- unconstrain_pars(sf3, get_inits(sf3)[[1]])
  checkEquals(i_u, rep(0, 2))
} 

.tearDown <- function() { 
  csv_fname <- 'teststanfit*.csv'
  system(paste('rm -rf ', csv_fname))
} 
