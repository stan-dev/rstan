print.conv.diag <-
function(x, ...){
  cat("Geweke Diag. & Heidelberger and Welch's Diag.\n")
  #print(x$sol_geweke)
  #cat("Heidelberger and Welch's Diag.\n")
  #print(x$sol_heidel)
  print(format(cbind(x$sol_geweke, x$sol_heidel), nsmall = 4))
  if(x$pass_ind){
    cat('\n')
    cat("The Chain has converged.\n")
  }else{
    warning("The Chain may not have converged. Consider a longer burn-in/warmup, speeding up the convergence by setting conv_speedup = T, or modification of the model.\n", call. = F, immediate. = T)
  }
}
