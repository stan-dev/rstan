ssquares <-
function (y, design_matrix, assign, factor_index){
  n_f <- length(factor_index)
  factor_SS <- array(NA, dim = c(1, n_f))
  # get the new design matrix according to assign and factor_index
  newassign <- array(0) #include intercept
  index_of_factors <- array(1) #include intercept
  for (i in 1:n_f){
    newassign <- c(newassign, assign[assign == factor_index[i]])
    index_of_factors <- c(index_of_factors, which(assign == factor_index[i]))
  }
  newdesign_matrix <- design_matrix[, index_of_factors]
  if (length(index_of_factors) > 1) 
    if (qr(newdesign_matrix)$rank < length(index_of_factors)) stop ('Colinearity in level 2 design matrix, model is unidentified!')
  #SS_TO <- sum((y - mean(y))^2) # total sum of squares for vectors
  centered.y <- scale(y, scale = F)
  SS_TO <- colSums(centered.y^2)
  full <- lsfit(newdesign_matrix, y, intercept = F)
  SSE <- colSums((full$residuals)^2)
  SS_full <- SS_TO - SSE # sum of squares of the full model
  if (n_f == 1){
    factor_SS[1] <- round(mean(SS_full), digits = 5)
      #paste(round(SS_full, digits = 5), ' (', round((SS_full)/SS_TO*100, digits = 2), '%)', sep="")
  }else{
    for (i in 1:n_f){
      factor_index_i <- which(newassign == factor_index[i])
      design_matrix_i <- newdesign_matrix[, -factor_index_i]
      model_i <- lsfit(design_matrix_i, y, intercept = F)
      SS_i <- SS_TO - colSums((model_i$residuals)^2)
      factor_SS[i] <- round(mean(SS_full - SS_i), digits = 5)
        #paste(round(SS_full - SS_i, digits = 5), ' (',round((SS_full - SS_i)/SS_TO*100, digits = 2),'%)',sep = '')
    }
  }
  return(list(factor_SS = factor_SS, SS_TO = round(mean(SS_TO), digits = 5), SSE = round(mean(SSE), digits = 5)))
}
