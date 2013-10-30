# test some functionality of function stan()
test_stan_fun_args <- function() {
  csv_fname <- 'tsfa.csv'
  csv_fname2 <- 'tsfa2.csv'
  diag_fname2 <- 'diag2.csv'
  model_code <- "
    parameters { 
      real y[2];
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
  "
  fit <- stan(model_code = model_code, 
              iter = 100, chains = 1, thin = 3, 
              sample_file = csv_fname)
  fit2 <- stan(fit = fit, iter = 1001, chains = 3, thin = 2,
               sample_file = csv_fname2, diagnostic_file = diag_fname2)
  fit3 <- stan(fit = fit, iter = 1001, chains = 3, thin = 2,
               control = list())
  fit4 <- stan(fit = fit, iter = 1001, chains = 3, thin = 2,
               control = list(adapt_gamma = .7)) 
} 

.tearDown <- function() { 
  csv_fname <- 'tsfa*.csv'
  csv_fname2 <- 'tsfa2*.csv'
  diag_fname2 <- 'diag2*.csv'
  system(paste('rm -rf ', csv_fname))
  system(paste('rm -rf ', csv_fname2))
  system(paste('rm -rf ', diag_fname2))
} 
