print.BANOVA.ordMultinomial <-
function(x, ...){
  cat('Call:\n')
  print(x$call)
  cat('\n Coefficients: \n')
  #print(x$coef.tables$full_table)
  #print(x$coef.tables$cutp_table)               
  print(data.frame(rbind(x$coef.tables$full_table, x$coef.tables$cutp_table)))
  
}
