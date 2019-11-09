print.BANOVA.Bernoulli <-
function(x, ...){
  cat('Call:\n')
  print(x$call)
  cat('\n Coefficients: \n')
  print(data.frame(x$coef.tables$full_table))
  
}
