traceplot <-
function (samples_l2_param, X_names = " ", Z_names, save = FALSE){
  row_names <- array(NA, dim = c(ncol(samples_l2_param),1))
  for (i in 1:length(X_names))
    for (j in 1:length(Z_names)){
      temp <- (i - 1)*length(Z_names) + j
      if (X_names[i] == " ")
        row_names[temp] <- Z_names[j]
      else
        row_names[temp] <- paste(X_names[i]," : ", Z_names[j])
    }
  plot_data <- samples_l2_param
  colnames(plot_data) <- row_names
  if (save == TRUE){
    pdf('traceplot.pdf')
    plot(mcmc(plot_data))
    dev.off()
  }else
    plot(mcmc(plot_data))
}
