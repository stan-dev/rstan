# to test part of stan_args.hpp
.setUp <- function() {
  require(rstan) 
  src <- ' 
    BEGIN_RCPP
    Rcpp::List lst(x); 
    rstan::stan_args args(lst); 
    return args.stan_args_to_rlist(); 
    END_RCPP
  ' 
  fx <- cxxfunction(signature(x = "list"), 
                    body = src, 
                    includes = "#include <rstan/stan_args.hpp>", 
                    plugin = "rstan", verbose = TRUE)
  assign("fx", fx, envir = .GlobalEnv)
}

test_stan_args_hppb <- function() {
  b1 <- fx(list(iter = 100, seed = 12354, method = 'optim')) 
  checkEquals(b1$random_seed, "12354")
}

test_stan_args_hpp <- function() { 

  a1 <- fx(list(iter = 100, thin = 100)) 
  a2 <- fx(list(iter = 100, thin = 3)) 
  a3 <- fx(list(iter = 5, thin = 3, refresh = -1, seed = "12345")) 
  a4 <- fx(list(iter = 5, thin = 3, refresh = -1, method = 'test_grad'))
  a4b <- fx(list(iter = 5, thin = 3, refresh = -1, method = 'test_grad',
                 control = list(epsilon = 1.3, error = 0.1)))
  a4c <- fx(list(iter = 5, thin = 3, refresh = -1, method = 'test_grad',
                 control = list(epsilon = 1.3)))
  a4d <- fx(list(iter = 5, thin = 3, refresh = -1, method = 'test_grad',
                 control = list(error = 1.3)))
  a5 <- fx(list(iter = 5, thin = 3, algorithm = 'HMC', 
                control = list(stepsize = .1)))
  a6 <- fx(list(iter = 5, thin = 3, algorithm = 'NUTS', 
                control = list(stepsize = .1, metric = 'unit_e')))
  a7 <- fx(list(iter = 5, thin = 3, algorithm = 'NUTS', 
                control = list(stepsize = .1, metric = 'diag_e')))
  a7b <- fx(list(iter = 5, thin = 3, algorithm = 'NUTS', 
                 control = list(stepsize = .1, metric = 'diag_e', 
                                adapt_term_buffer = 4, 
                                adapt_window = 30,
                                adapt_init_buffer = 40)))
  a8 <- fx(list(iter = 5, thin = 3, algorithm = 'NUTS', 
                control = list(stepsize = .1, metric = 'dense_e')))
  a9 <- fx(list(iter = 5, thin = 3, algorithm = 'HMC', 
                control = list(stepsize = .1, metric = 'unit_e')))
  a10 <- fx(list(iter = 5, thin = 3, algorithm = 'HMC', 
                 control = list(stepsize = .1, metric = 'diag_e')))
  a11 <- fx(list(iter = 5, thin = 3, algorithm = 'HMC', 
                 control = list(stepsize = .1, metric = 'dense_e')))

  checkEquals(a1$iter, 100) 
  checkEquals(a1$thin, 100) 
  checkEquals(a2$iter, 100) 
  checkEquals(a2$thin, 3) 
  checkEquals(a2$control$adapt_init_buffer, 75)
  checkEquals(a2$control$adapt_term_buffer, 50)
  checkEquals(a2$control$adapt_window, 25)
  checkEquals(a3$iter, 5) 
  checkEquals(a3$random_seed, "12345")
  checkEquals(a3$test_grad, FALSE)
  checkEquals(a3$refresh, -1) 
  checkEquals(a3$control$adapt_init_buffer, 75)
  checkEquals(a3$control$adapt_term_buffer, 50)
  checkEquals(a3$control$adapt_window, 25)
  checkEquals(a4$test_grad, TRUE)
  checkEquals(a4$method, "test_grad")
  checkEquals(a4$control$epsilon, 1e-6)
  checkEquals(a4$control$error, 1e-6)
  checkEquals(a4b$control$epsilon, 1.3)
  checkEquals(a4b$control$error, 1e-1)
  checkEquals(a4c$control$epsilon, 1.3)
  checkEquals(a4c$control$error, 1e-6)
  checkEquals(a4c$control$error, 1e-6)
  checkEquals(a4c$control$epsilon, 1.3)
  # checkTrue()
  checkException(fx(list(iter = 5, thin = 3, refresh = -1, test_grad = TRUE, seed = '111111111111111111111'))) 
  checkEquals(a6$control$metric, "unit_e")
  checkEquals(a7$control$metric, "diag_e")
  checkEquals(a7b$control$adapt_init_buffer, 40)
  checkEquals(a7b$control$adapt_term_buffer, 4)
  checkEquals(a7b$control$adapt_window, 30)
  checkEquals(a8$control$metric, "dense_e")
  checkEquals(a9$control$metric, "unit_e")
  checkEquals(a10$control$metric, "diag_e")
  checkEquals(a11$control$metric, "dense_e")
  checkEquals(a11$control$adapt_window, 25)
} 

