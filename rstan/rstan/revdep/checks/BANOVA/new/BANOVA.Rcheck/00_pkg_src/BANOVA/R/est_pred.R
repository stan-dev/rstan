est_pred <-
function(f, est_matrix, l1_matrix, l2_matrix, cutp0, cutp1){
  temp <- l1_matrix %*% est_matrix %*% t(l2_matrix)
  if (sum(cutp0 != cutp1) == 0){ 
    est <- 1 - f(temp)
  }else if(sum(cutp1 != 0) == 0){
    temp <- temp - cutp0
    est <- f(temp) 
  }else{
    est_cutp0 <- temp - cutp0
    est_cutp1 <- temp - cutp1
    est <- f(est_cutp0) - f(est_cutp1)
  }
  return(est)
}
