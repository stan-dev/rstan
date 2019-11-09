multi.predict.means <-
function(samples_l2_param, dataX, dataZ, X, X_original_choice, Z_full, mf1, mf2, Xsamples = NULL, Zsamples = NULL){
  if(is.null(Xsamples)){ 
    Xsamples <- dataX
  }  
  if(is.null(Zsamples)){ 
    Zsamples <- dataZ
  }
  if(!is.null(mf2))
    if (length(Xsamples) != nrow(Zsamples)) stop('X, Z samples dimension mismatch.')
  if (class(Xsamples) != 'list') stop('Xsamples must be a list of features data.')
  for (i in 1:length(Xsamples))
    if (ncol(Xsamples[[i]]) != ncol(dataX[[1]]) || nrow(Xsamples[[i]]) != nrow(dataX[[1]]))
      stop('level 1 samples dimension mismatch!')
  #if (is.vector(samples)) samples <- matrix(samples,  nrow  = 1)
  #if (is.vector(X)) X <- matrix(X,  ncol = 1)
  if (is.vector(Z_full)) Z_full <- matrix(Z_full,  ncol = 1)
  if(!is.null(mf2))
    if (ncol(Zsamples) != ncol(dataZ)) stop('level 2 samples dimension mismatch!')
  
  n_iter <- nrow(samples_l2_param)
  num_l1 <- ncol(X[[1]])
  if(is.null(mf2)){
    num_l2 <- 1
  }else{
    num_l2 <- ncol(Z_full)
  }
  est_matrix <- array(0 , dim = c(num_l1, num_l2, n_iter))
  for (i in 1:num_l1){
    for (j in 1:n_iter)
      est_matrix[i,,j] <- samples_l2_param[j,((i-1)*num_l2+1):((i-1)*num_l2+num_l2)]
  }
  
  l1_names <- attr(mf1, 'names')
  if(is.null(mf2)){
    l2_names <- c(" ")
  }else{
    l2_names <- attr(mf2, 'names')
  }
  l1_index_in_data <- which(colnames(dataX[[1]]) %in% l1_names)
  l2_index_in_data <- which(colnames(dataZ) %in% l2_names)
  # find index of level 1 factors and numeric variables
  l1_factor_index_in_data <- array(dim = 0)
  l1_numeric_index_in_data <- array(dim = 0)
  if (length(l1_index_in_data) > 0){
    for (i in 1: length(l1_index_in_data)){
      if (class(dataX[[1]][,l1_index_in_data[i]]) == 'factor')
        l1_factor_index_in_data <- c(l1_factor_index_in_data, l1_index_in_data[i])
      if (class(dataX[[1]][,l1_index_in_data[i]]) == 'integer' || class(dataX[[1]][,l1_index_in_data[i]]) == 'numeric')
        l1_numeric_index_in_data <- c(l1_numeric_index_in_data, l1_index_in_data[i])
    }
  }
  # find index of level 2 factors and numeric variables
  l2_factor_index_in_data <- array(dim = 0)
  l2_numeric_index_in_data <- array(dim = 0)
  if (length(l2_index_in_data) > 0){
    for (i in 1: length(l2_index_in_data)){
      if (class(dataZ[,l2_index_in_data[i]]) == 'factor')
        l2_factor_index_in_data <- c(l2_factor_index_in_data, l2_index_in_data[i])
      if (class(dataZ[,l2_index_in_data[i]]) == 'integer' || class(dataZ[,l2_index_in_data[i]]) == 'numeric')
        l2_numeric_index_in_data <- c(l2_numeric_index_in_data, l2_index_in_data[i])
    }
  }
  
  n_choice <- length(X)
  prediction <- matrix(0, nrow = n_choice*length(Xsamples), ncol = 5)
  colnames(prediction) <- c('Sample number', 'Choice','Median', '2.5%', '97.5%')
  est_samples <- matrix(0, nrow = n_choice, ncol = n_iter)
  temp_index <- 1
  for (n_sample in 1 : length(Xsamples)){
    l1_vector <- list()
    for (n_c in 1: n_choice){
      l1_vector[[n_c]] <- matrix(rep(0, num_l1), nrow = 1)
      if (n_c != 1) l1_vector[[n_c]][1,n_c-1] <- 1
      if (length(l1_factor_index_in_data) > 0){
        index_row <- rowMatch(Xsamples[[n_sample]][n_c, l1_factor_index_in_data], X_original_choice[[n_c]][, l1_factor_index_in_data])    
        if (is.na(index_row)) stop('Bad samples provided! Could not find matching factors!')
        l1_vector[[n_c]] <- X[[n_c]][index_row, ]
        
      }
      if (length(l1_numeric_index_in_data) > 0)
        for (i in 1:length(l1_numeric_index_in_data)){
          l1_vector[[n_c]][attr(X[[1]], 'numeric_index')[-c(1:(n_choice - 1))][i]] <- Xsamples[[n_sample]][n_c, l1_numeric_index_in_data[i]]
        }
    }
    l2_vector <- matrix(c(1, rep(0, num_l2-1)), nrow = 1)
    if (length(l2_factor_index_in_data) > 0){
      index_row_Z <- rowMatch(Zsamples[n_sample, l2_factor_index_in_data], dataZ[, l2_factor_index_in_data])
      if (is.na(index_row_Z)) stop('Bad samples provided! Could not find matching factors!')
      l2_vector <- Z_full[index_row_Z, ]
    }
    
    if (length(l2_numeric_index_in_data) > 0)
      for (i in 1:length(l2_numeric_index_in_data))
        l2_vector[attr(Z_full, 'numeric_index')[i]] <- Zsamples[n_sample,l2_numeric_index_in_data[i]]
    for(n_c in 1: n_choice){
      for (n_i in 1:n_iter){
        if (class(est_matrix[,,n_i]) == 'numeric' | class(est_matrix[,,n_i]) == 'integer'){ # not a matrix, R somehow automatically convert dim(1,n) matrix to a vector
          if (length(l1_vector[[n_c]]) == 1) temp <- matrix(est_matrix[,,n_i], nrow = 1)
          if (length(l2_vector) == 1) temp <- matrix(est_matrix[,,n_i], ncol = 1)
          #print(l1_vector[[n_c]])
          est_samples[n_c, n_i] <- exp(matrix(l1_vector[[n_c]], nrow = 1) %*% temp %*% t(matrix(l2_vector, nrow = 1)))
          
        }else{
          est_samples[n_c, n_i] <- exp(matrix(l1_vector[[n_c]], nrow = 1) %*% est_matrix[,,n_i] %*% t(matrix(l2_vector, nrow = 1)))
        }
      }
    }
    est_prob_sum <- matrix(0, nrow = 1, ncol = n_iter)
    for (n_c in 1:n_choice)
      est_prob_sum <- est_prob_sum + est_samples[n_c,]
    est_prob <- matrix(0, nrow = n_choice, ncol = n_iter)
    for (n_c in 1:n_choice){
      est_prob[n_c, ] <- est_samples[n_c, ] / est_prob_sum
      means <- median(est_prob[n_c, ])
      quantile_025975 <- quantile(est_prob[n_c, ], probs = c(0.025, 0.975))
      prediction[temp_index, ] <- c(n_sample, n_c, round(c(means, quantile_025975), digits = 5))
      temp_index <- temp_index + 1
    }
  }
  return(prediction)
}
