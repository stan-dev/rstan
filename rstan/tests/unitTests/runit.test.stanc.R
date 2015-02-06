test_long_erromsg <- function() {
  code <- "parameters { cov_matrix[3] y; }  
           model { y ~ normal(0, 1); }"
  checkException(r <- stanc(model_code = code))
  warning.length <- getOption("warning.length")
  emsg <- geterrmessage()
  checkTrue(nchar(emsg) > 1000)
  checkEquals(warning.length, getOption("warning.length"))
} 
