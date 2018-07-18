test_getting_diagnostics <- function() {
  exfit <- read_stan_csv("unitTests/test_fit_diagnostics.csv")
  checkEqualsNumeric(get_bfmi(exfit), 1.17558932146879)
  checkEquals(get_low_bfmi_chains(exfit), integer(0))
  
  expected_divergent <- logical(10)
  expected_divergent[c(2,3)] <- TRUE
  checkEquals(get_divergent_iterations(exfit), expected_divergent)
  checkEquals(get_num_divergent(exfit), 2)
  
  expected_max_treedepth <- logical(10)
  expected_max_treedepth[c(1,5)] <- TRUE
  checkEquals(get_max_treedepth_iterations(exfit), expected_max_treedepth)
  checkEquals(get_num_max_treedepth(exfit), 2)
  
  checkEquals(get_num_leapfrog_per_iteration(exfit), c(1,3,5,5,7,5,3,1,1,7))  
}
