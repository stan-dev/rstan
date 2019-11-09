table.coefficients <-
function (samples_l2_param, l2_parameters, X_names, Z_names, X_assign, Z_assign, samples_cutp_param = array(dim = 0)){
  n_p <- length(l2_parameters)
  row_names <- array(NA, dim = c(n_p,1))
  row_indices <- matrix(NA, nrow = length(X_names)*length(Z_names), ncol = 2) # store indices (corresponding to level 1 and 2 variables) of each row in coefficient table, used for p-value table
  for (i in 1:length(X_names))
    for (j in 1:length(Z_names)){
      temp <- (i - 1)*length(Z_names) + j
      if (X_names[i] == " "){
        row_names[temp] <- Z_names[j]
      }else{
        row_names[temp] <- paste(X_names[i]," : ", Z_names[j])
      }
      row_indices[temp, 1] <- X_assign[i]
      row_indices[temp, 2] <- Z_assign[j]
    }
  result_table <- data.frame(matrix(NA, nrow = n_p, ncol = 1+1+2+1+1)) # the table stores statistics of coefficients 
  rownames(result_table) <- row_names
  colnames(result_table) <- c("mean","SD","Quantile0.025","Quantile0.975","p-value","Signif.codes")
  for (i in 1:n_p){
    result_table[i,1] <- round(mean(samples_l2_param[,i]), digits = 4)
    result_table[i,2] <- round(as.numeric(sd(samples_l2_param[,i])), digits = 4)
    result_table[i,3:4] <- round(as.numeric(quantile(samples_l2_param[,i],c(0.025,0.975))), digits = 4)
  }
  p_values <- pValues(samples_l2_param)
  #result_table <- round(result_table, digits = 4)
  result_table[,5] <- round(p_values, digits = 4)
  coeff_table <- result_table[, c(1,3,4,5)]
  #result_table <- round(result_table, digits = 4)
  result_table <- format(result_table, nsmall = 4)
  for (i in 1:n_p){
    if (p_values[i] <= 0.001) {
      result_table[i,6] <- '***'
      if (p_values[i] == 0) result_table[i,5] <- '<0.0001'
    }
    else if (p_values[i] <= 0.01 & p_values[i] > 0.001 ) result_table[i,6] <- '**'
    else if	(p_values[i] <= 0.05 & p_values[i] > 0.01 ) result_table[i,6] <- '*'
    else if	(p_values[i] <= 0.1 & p_values[i] > 0.05 ) result_table[i,6] <- '.'
    else if	(p_values[i] <= 1 & p_values[i] > 0.1 ) result_table[i,6] <- ' '
  }

  # for multinomial cutpoint model
  if (length(samples_cutp_param) != 0){
    rn_cutpresults <- array(NA,ncol(samples_cutp_param))
    for (i in 1:ncol(samples_cutp_param))
      rn_cutpresults[i] <-  paste('Cutpoint[',i+1,']', sep="")
    cutpresults <- matrix(NA,nrow = ncol(samples_cutp_param), ncol = 6)
    rownames(cutpresults) <- rn_cutpresults
    colnames(cutpresults) <- colnames(result_table)
    for (i in 1:ncol(samples_cutp_param)){
      cutpresults[i,1] <- round(mean(samples_cutp_param[,i]),4)
      cutpresults[i,2] <- round(sd(samples_cutp_param[,i]),4)
      cutpresults[i,3:4] <- round(quantile(samples_cutp_param[,i],c(0.025,0.975)),4)
    }
  }
  #full_table <- as.table(result_table)
  results <- list()
  results$row_indices <- row_indices
  results$full_table <- result_table
  results$coeff_table <- data.frame(coeff_table)
  if (length(samples_cutp_param) != 0)
    results$cutp_table <- cutpresults
  return(results) # return the coefficient table with mean and quantiles for table of means computation
}
