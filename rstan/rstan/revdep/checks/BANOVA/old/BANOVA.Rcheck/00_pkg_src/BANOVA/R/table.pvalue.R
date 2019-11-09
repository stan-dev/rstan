table.pvalue <-
function (coeff_table, row_indices, l1_names, l2_names){
  if (length(l1_names) == 1 && length(l2_names) == 1) coeff_table <- data.frame(matrix(coeff_table, nrow = 1)) # the case that there is only one intercept 
  n_l1_p <- length(l1_names)
  n_l2_p <- length(l2_names)
  ptable <- data.frame(matrix(NA, nrow =  n_l1_p, ncol = n_l2_p))
  for (i in 1:n_l1_p)
    for (j in 1:n_l2_p){
      # find rows in coeff_table corresponding to each factor or interaction
      rows <- which(row_indices[,1] == i & row_indices[,2] == j)
      # find the min and max of the means among these rows
      temp_mean <- coeff_table[rows, 1]
      temp_pvalues <- coeff_table[rows, 4]
      # find the p-values of corresponding variables with min and max mean 
      p_min <- temp_pvalues[which.min(temp_mean)]
      p_max <- temp_pvalues[which.max(temp_mean)]
      ptable[i, j] <- round(min(p_min, p_max), digits = 4)
      if (ptable[i, j] == 0) ptable[i, j] <- '<0.0001' else ptable[i, j] <- format(as.numeric(ptable[i, j]), nsmall = 4)
    }
  rownames(ptable) <- l1_names
  colnames(ptable) <- l2_names
  return(ptable)
}
