print.summary.BANOVA <- function(x, ...){
  cat('Call:\n')
  print(x$call)
  
  cat('\nConvergence diagnostics:\n')
  print(x$conv)
  
  cat('\nTable of sum of squares & effect sizes:\n')
  if (length(x$anova.table) >2){
    for (i in 1:length(x$anova.table)){
      cat('\nChoice: ', i, '\n')
      print(x$anova.table[[i]])
    }
  }else{
    print(x$anova.table)
  }
  cat('\nTable of p-values (Multidimensional): \n')
  print(x$pvalue.table)
  
  cat('\nTable of coefficients: \n')
  printCoefmat(x$coef.table)
  
  if (!is.null(x$R2)){
    cat('\nMultiple R-squared: ', x$R2, '\n')
  }
    
  if (!((x$model_name == "BANOVA.Multinomial")&&(x$single_level))){
    cat('\nTable of predictions: \n')
    table.predictions(x$full_object)
  }
}