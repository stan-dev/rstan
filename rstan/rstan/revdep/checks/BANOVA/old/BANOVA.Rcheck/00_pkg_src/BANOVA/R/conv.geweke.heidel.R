conv.geweke.heidel <-
function (samples_l2_param, X_names, Z_names){
  row_names <- array(NA, dim = c(ncol(samples_l2_param),1))
  for (i in 1:length(X_names))
    for (j in 1:length(Z_names)){
      temp <- (i - 1)*length(Z_names) + j
      if (X_names[i] == " ")
        row_names[temp] <- Z_names[j]
      else
        row_names[temp] <- paste(X_names[i]," : ", Z_names[j])
    }
  pass_ind <- T
  sol_converge <- geweke.diag(as.mcmc(samples_l2_param),0.2,0.5)
  p_values_converge <- 2*pnorm(-abs(sol_converge$z)) 
  decision_conerge <- array(NA,length(p_values_converge)) 
  decision_conerge[p_values_converge>=0.001] <- 'passed' 
  decision_conerge[p_values_converge<0.001] <- 'failed' 
  if (sum(p_values_converge<0.001) > 0) pass_ind <- F
  sol_geweke <- data.frame('Geweke stationarity test' = decision_conerge, 'Geweke convergence p value' = round(p_values_converge, digits = 4))
  #sol_geweke <- data.frame(cbind(decision_conerge, as.numeric(round(p_values_converge, digits = 4))))
  #print(str(sol_geweke))
  rownames(sol_geweke) <- row_names
  colnames(sol_geweke) <- c('Geweke stationarity test','Geweke convergence p value') 
  
  hddata <- samples_l2_param
  colnames(hddata) <- row_names
  sol_heidel_temp <- heidel.diag(as.mcmc(hddata))
  sol_heidel <- data.frame(H_W_stat = sol_heidel_temp[, 1], H_W_p = round(sol_heidel_temp[, 3], 4))
  colnames(sol_heidel) <- c('H. & W. stationarity test', 'H. & W. convergence p value')
  #print(sol_heidel)
  #sol_heidel[,3] <- round(sol_heidel[,3],4)
  #sol_heidel[,3] <- round(sol_heidel[,3],4)
  #sol_heidel[,5] <- round(sol_heidel[,5],4)
  #sol_heidel[,6] <- round(sol_heidel[,6],4)

  if (sum(sol_heidel[,1]<1) > 0) pass_ind <- F
  sol_heidel[which(sol_heidel[,1]==1),1]<-'passed'
  sol_heidel[which(sol_heidel[,1]<1),1]<-'failed'

  sol <- list(sol_geweke = sol_geweke, sol_heidel = sol_heidel, pass_ind = pass_ind)
}
