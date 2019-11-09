BANOVA.T <-
function(l1_formula = 'NA', l2_formula = 'NA', data, id, l1_hyper = c(1, 1, 1),l2_hyper = c(1,1,0.0001), burnin = 5000, sample = 2000, thin = 10, adapt = 0, conv_speedup = F, jags = runjags.getOption('jagspath')){
  # if (jags == 'JAGS not found') stop('Please install the JAGS software (Version 3.4 or above). Check http://mcmc-jags.sourceforge.net/')
  sol <- BANOVA.TNormal(l1_formula, l2_formula, data, id, l1_hyper = l1_hyper, l2_hyper = l2_hyper, burnin = burnin, sample = sample, thin = thin, adapt = adapt, conv_speedup = conv_speedup, jags = jags)
  sol$call <- match.call()
  class(sol) <- 'BANOVA.T'
  sol
}
