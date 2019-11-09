table.ANOVA <-
function(samples_l1_param, X, Z){
  if (length(attr(Z, 'varNames')) > 1){
    num_l1_v <- ncol(X)
    #make sure to exclude the numeric variables
    numeric_index <- which(attr(Z,'dataClasses') == 'numeric' | attr(Z,'dataClasses') == 'integer')
    factor_index <- which(attr(Z,'dataClasses') == 'factor')
    temp_names <- attr(attr(Z,'dataClasses'), 'names')
    numeric_names <- temp_names[numeric_index]
    numeric_index_in_Vnames <- array(dim = 0)
    if (length(numeric_index) > 0)
      for(i in 1:length(numeric_index))
        numeric_index_in_Vnames <- c(numeric_index_in_Vnames, which(attr(Z, 'varNames') == numeric_names[i]))
    #num_l2_v <- length(attr(Z, 'varNames')) - 1 - length(numeric_index) # exclude the intercept and numeric variables
    num_l2_v <- length(factor_index) + length(attr(Z,'interactions_index'))
    if (num_l2_v > 0){
      num_id <- nrow(Z)
      anova_table <- matrix(NA, nrow = num_l1_v, ncol = num_l2_v+2) 
      rownames(anova_table) <- colnames(X)
      #colnames(anova_table) <- c(attr(Z, 'varNames')[-c(1, numeric_index_in_Vnames)], 'Residuals', 'Total') # exclude the intercept and numerics
      temp_full_names <- attr(Z, 'varNames')[-1] # exclude intercept
      colnames(anova_table) <- c(temp_names[factor_index], temp_full_names[attr(Z,'interactions_index')], 'Residuals', 'Total')
      assign <- attr(Z, 'assign')
      Z_factor_index <- c(factor_index, attr(Z,'interactions_index'))
      for (i in 1:num_l1_v){
        y <- array(0, dim = c(num_id, 1))
        for (j in 1:num_id)
          y[j] <- mean(samples_l1_param[,(i-1)*num_id + j])
        SS <- ssquares(y, Z, assign, Z_factor_index)
        anova_table[i,] <- c(SS$factor_SS, SS$SSE, SS$SS_TO)
      }
      return(data.frame(anova_table))
    }else
      return(NA)
  }else
    return(NA)
}
