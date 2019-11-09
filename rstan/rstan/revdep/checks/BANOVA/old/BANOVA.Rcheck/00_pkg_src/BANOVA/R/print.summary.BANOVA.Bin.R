print.summary.BANOVA.Binomial <-
function(x, ...){
  cat('Call:\n')
  print(x$call)
  
  cat('\nConvergence diagnostics:\n')
  print(x$conv)
  
  print(x$anova.table)
  
  cat('\nTable of p-values (Multidimensional): \n')
  print(x$pvalue.table)
  
  cat('\nTable of coefficients: \n')
  printCoefmat(x$coef.table)
  
  cat('\nTable of predictions: \n')
  table.predictions(x$full_object)
  
}
