pValues <-
function (sample){
  n_row <- nrow(sample)
  n_col <- ncol(sample)
  p_value <- array(NA,n_col)
  for (i in 1:n_col){
    n_tailp <- length(which(sample[,i] > 0))
    n_tailn <- length(which(sample[,i] < 0))
    p_value[i] <- round(min(n_tailp,n_tailn) / n_row, digits = 5)

  }
  return (2 * p_value) # two tailed test
}
