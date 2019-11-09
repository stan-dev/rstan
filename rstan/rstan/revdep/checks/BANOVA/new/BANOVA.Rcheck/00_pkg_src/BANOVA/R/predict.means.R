predict.means <-
function (samples_l2_param, data, X, Z_full, mf1, mf2 = NULL, samples = NULL, samples_cutp_param = NA, model = NA, l2_sd = NA, n_trials = NULL){
  if(is.null(samples)) samples <- data
  if (is.vector(samples)) samples <- matrix(samples,  nrow  = 1)
  if (is.vector(X)) X <- matrix(X,  ncol = 1)
  if (is.vector(Z_full)) X <- matrix(Z_full,  ncol = 1)
  if (ncol(samples) != ncol(data)) stop("Samples' dimension mismatch!")
  
  n_iter <- nrow(samples_l2_param)
  num_l1 <- ncol(X)
  if(is.null(mf2) || is.null(Z_full)){
    num_l2 <- 1
  }else{
    num_l2 <- ncol(Z_full)
  }
  
  est_matrix <- array(0 , dim = c(num_l1, num_l2, n_iter))
  for (i in 1:num_l1){
    for (j in 1:n_iter)
      est_matrix[i,,j] <- samples_l2_param[j,((i-1)*num_l2+1):((i-1)*num_l2+num_l2)]
  }
  
  if (model == 'NormalNormal'){
    link_inv <- identity
  }else if (model == 'PoissonNormal'){
    link_inv <- exp
    l2_var <- l2_sd^2
  }else if (model == 'BernNormal' || model == 'MultinomialordNormal'){
    link_inv <- function(x) return(exp(x)/(exp(x) + 1))
  }
  
  l1_names <- attr(mf1, 'names')[-1] # exclude y
  if (is.null(mf2)){
    l2_names <- c(" ")
  }else{
    l2_names <- attr(mf2, 'names')
  }
  l1_index_in_data <- which(colnames(data) %in% l1_names)
  l2_index_in_data <- which(colnames(data) %in% l2_names)
  # find index of level 1 factors and numeric variables
  l1_factor_index_in_data <- array(dim = 0)
  l1_numeric_index_in_data <- array(dim = 0)
  if (length(l1_index_in_data) > 0){
    for (i in 1: length(l1_index_in_data)){
      if (class(data[,l1_index_in_data[i]]) == 'factor')
        l1_factor_index_in_data <- c(l1_factor_index_in_data, l1_index_in_data[i])
      if (class(data[,l1_index_in_data[i]]) == 'integer' || class(data[,l1_index_in_data[i]]) == 'numeric')
        l1_numeric_index_in_data <- c(l1_numeric_index_in_data, l1_index_in_data[i])
    }
  }
  # find index of level 2 factors and numeric variables
  l2_factor_index_in_data <- array(dim = 0)
  l2_numeric_index_in_data <- array(dim = 0)
  if (length(l2_index_in_data) > 0){
    for (i in 1: length(l2_index_in_data)){
      if (class(data[,l2_index_in_data[i]]) == 'factor')
        l2_factor_index_in_data <- c(l2_factor_index_in_data, l2_index_in_data[i])
      if (class(data[,l2_index_in_data[i]]) == 'integer' || class(data[,l2_index_in_data[i]]) == 'numeric')
        l2_numeric_index_in_data <- c(l2_numeric_index_in_data, l2_index_in_data[i])
    }
  }
  l12_factor_index <- c(l1_factor_index_in_data, l2_factor_index_in_data) 
  if (model == 'BernNormal'){
    y_pred <- matrix(0, nrow = nrow(samples), ncol = 3)
    colnames(y_pred) <- c('Median', '2.5%', '97.5%')
    est_samples <- matrix(0, nrow = 1, ncol = n_iter)
    sample_size = mean(n_trials)
    for (n_sample in 1 : nrow(samples)){
      l1_vector <- matrix(c(1, rep(0, num_l1-1)), nrow = 1)
      l2_vector <- matrix(c(1, rep(0, num_l2-1)), nrow = 1)
      if (length(l12_factor_index) > 0){
        index_row <- rowMatch(samples[n_sample, l12_factor_index], data[, l12_factor_index])    
        if (is.na(index_row)) stop('Bad samples provided! Could not find matching factors!')
        l1_vector <- X[index_row, ]
        if (!is.null(Z_full))
          l2_vector <- Z_full[index_row, ]
      }
      # TODO numeric variables included in prediction +-SD
      #if (length(l1_numeric_index_in_data) > 0)
      #  for (i in 1:length(l1_numeric_index_in_data))
      #    l1_vector[attr(X, 'numeric_index')[i]] <- samples[n_sample,l1_numeric_index_in_data[i]]
      
      #if (length(l2_numeric_index_in_data) > 0)
      #  for (i in 1:length(l2_numeric_index_in_data))
      #    l2_vector[attr(Z_full, 'numeric_index')[i]] <- samples[n_sample,l2_numeric_index_in_data[i]]
      if (sum(is.na(l2_sd)) > 0){
        for (n_i in 1:n_iter){
          if (class(est_matrix[,,n_i]) == 'numeric' | class(est_matrix[,,n_i]) == 'integer'){ # not a matrix, R somehow automatically convert dim(1,n) matrix to a vector
            if (length(l1_vector) == 1) temp <- matrix(est_matrix[,,n_i], nrow = 1)
            if (length(l2_vector) == 1) temp <- matrix(est_matrix[,,n_i], ncol = 1)
            est_samples[n_i] <- matrix(l1_vector, nrow = 1) %*% temp %*% t(matrix(l2_vector, nrow = 1))
          }else{
            est_samples[n_i] <- matrix(l1_vector, nrow = 1) %*% est_matrix[,,n_i] %*% t(matrix(l2_vector, nrow = 1))
          }
        }
      }else{
        for (n_i in 1:n_iter){
          if (class(est_matrix[,,n_i]) == 'numeric' | class(est_matrix[,,n_i]) == 'integer'){ # not a matrix, R somehow automatically convert dim(1,n) matrix to a vector
            if (length(l1_vector) == 1) temp <- matrix(est_matrix[,,n_i], nrow = 1)
            if (length(l2_vector) == 1) temp <- matrix(est_matrix[,,n_i], ncol = 1)
            est_samples[n_i] <- matrix(l1_vector, nrow = 1) %*% temp %*% t(matrix(l2_vector, nrow = 1)) + sum(l2_var[n_i,])/2
          }else{
            est_samples[n_i] <- matrix(l1_vector, nrow = 1) %*% est_matrix[,,n_i] %*% t(matrix(l2_vector, nrow = 1)) + sum(l2_var[n_i,])/2
          }
        }
      }
      
      means <- apply(est_samples, 1, median) 
      quantile_025 <- apply(est_samples, 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
      quantile_975 <- apply(est_samples, 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
      y_pred[n_sample, 1:3] <- c(link_inv(means) * sample_size, link_inv(quantile_025) * sample_size, link_inv(quantile_975) * sample_size)
    }
    y_pred <- round(y_pred, digits = 4)
    return(y_pred)
  }
  if (model != 'MultinomialordNormal' && model != 'BernNormal'){
    y_pred <- matrix(0, nrow = nrow(samples), ncol = 3)
    colnames(y_pred) <- c('Median', '2.5%', '97.5%')
    est_samples <- matrix(0, nrow = 1, ncol = n_iter)
    for (n_sample in 1 : nrow(samples)){
      l1_vector <- matrix(c(1, rep(0, num_l1-1)), nrow = 1)
      l2_vector <- matrix(c(1, rep(0, num_l2-1)), nrow = 1)
      if (length(l12_factor_index) > 0){
        index_row <- rowMatch(samples[n_sample, l12_factor_index], data[, l12_factor_index])    
        if (is.na(index_row)) stop('Bad samples provided! Could not find matching factors!')
        l1_vector <- X[index_row, ]
        if (!is.null(Z_full))
          l2_vector <- Z_full[index_row, ]
      }
      # TODO numeric variables included in prediction
      #if (length(l1_numeric_index_in_data) > 0)
      #  for (i in 1:length(l1_numeric_index_in_data))
      #    l1_vector[attr(X, 'numeric_index')[i]] <- samples[n_sample,l1_numeric_index_in_data[i]]
      
      #if (length(l2_numeric_index_in_data) > 0)
      #  for (i in 1:length(l2_numeric_index_in_data))
      #    l2_vector[attr(Z_full, 'numeric_index')[i]] <- samples[n_sample,l2_numeric_index_in_data[i]]
      if (sum(is.na(l2_sd)) > 0){
        for (n_i in 1:n_iter){
          if (class(est_matrix[,,n_i]) == 'numeric' | class(est_matrix[,,n_i]) == 'integer'){ # not a matrix, R somehow automatically convert dim(1,n) matrix to a vector
            if (length(l1_vector) == 1) temp <- matrix(est_matrix[,,n_i], nrow = 1)
            if (length(l2_vector) == 1) temp <- matrix(est_matrix[,,n_i], ncol = 1)
            est_samples[n_i] <- matrix(l1_vector, nrow = 1) %*% temp %*% t(matrix(l2_vector, nrow = 1))
          }else{
            est_samples[n_i] <- matrix(l1_vector, nrow = 1) %*% est_matrix[,,n_i] %*% t(matrix(l2_vector, nrow = 1))
          }
        }
      }else{
        for (n_i in 1:n_iter){
          if (class(est_matrix[,,n_i]) == 'numeric' | class(est_matrix[,,n_i]) == 'integer'){ # not a matrix, R somehow automatically convert dim(1,n) matrix to a vector
            if (length(l1_vector) == 1) temp <- matrix(est_matrix[,,n_i], nrow = 1)
            if (length(l2_vector) == 1) temp <- matrix(est_matrix[,,n_i], ncol = 1)
            est_samples[n_i] <- matrix(l1_vector, nrow = 1) %*% temp %*% t(matrix(l2_vector, nrow = 1)) + sum(l2_var[n_i,])/2
          }else{
            est_samples[n_i] <- matrix(l1_vector, nrow = 1) %*% est_matrix[,,n_i] %*% t(matrix(l2_vector, nrow = 1)) + sum(l2_var[n_i,])/2
          }
        }
      }
      
      means <- apply(est_samples, 1, median)
      quantile_025 <- apply(est_samples, 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
      quantile_975 <- apply(est_samples, 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
      y_pred[n_sample, 1:3] <- c(link_inv(means), link_inv(quantile_025), link_inv(quantile_975))
    }
    y_pred <- round(y_pred, digits = 4)
    return(y_pred)
  }
  if (model == 'MultinomialordNormal'){
    n.cut <- ncol(samples_cutp_param) + 1
    cut_samples <- cbind(0, 0, samples_cutp_param, 0) # the first cut point is 0, the previous is set to be 0 to compute the prob. of Y == 1 (1 - logit^-1), the last 0 is for the prob. Y == K
    est_samples <- matrix(0, nrow = 1, ncol = n_iter)
    prediction <- matrix(0, nrow = (n.cut+1)*nrow(samples), ncol = 5)
    colnames(prediction) <- c('Sample number', 'Response','Median', '2.5%', '97.5%')
    temp_index <- 1
    for (n_sample in 1 : nrow(samples)){
      l1_vector <- matrix(c(1, rep(0, num_l1-1)), nrow = 1)
      l2_vector <- matrix(c(1, rep(0, num_l2-1)), nrow = 1)
      if (length(l12_factor_index) > 0){
        index_row <- rowMatch(samples[n_sample, l12_factor_index], data[, l12_factor_index])    
        if (is.na(index_row)) stop('Bad samples provided! Could not find matching factors!')
        l1_vector <- X[index_row, ]
        if (!is.null(Z_full))
          l2_vector <- Z_full[index_row, ]
      }
      
      # TO numeric variables included in prediction
      #if (length(l1_numeric_index_in_data) > 0)
      #  for (i in 1:length(l1_numeric_index_in_data))
      #    l1_vector[attr(X, 'numeric_index')[i]] <- samples[n_sample,l1_numeric_index_in_data[i]]
      
      #if (length(l2_numeric_index_in_data) > 0)
      #  for (i in 1:length(l2_numeric_index_in_data))
      #    l2_vector[attr(Z_full, 'numeric_index')[i]] <- samples[n_sample,l2_numeric_index_in_data[i]]
      
      for (y in 1:(n.cut + 1)){
        for (n_i in 1:n_iter){
          if (class(est_matrix[,,n_i]) == 'numeric' | class(est_matrix[,,n_i]) == 'integer'){ # not a matrix, R somehow automatically convert dim(1,n) matrix to a vector
            if (length(l1_vector) == 1) temp <- matrix(est_matrix[,,n_i], nrow = 1)
            if (length(l2_vector) == 1) temp <- matrix(est_matrix[,,n_i], ncol = 1)
            est_samples[n_i] <- est_pred(link_inv, temp, matrix(l1_vector, nrow = 1), matrix(l2_vector, nrow = 1), 
                                         cut_samples[n_i,y], cut_samples[n_i,y+1])
          }else{
            est_samples[n_i] <- est_pred(link_inv, est_matrix[,,n_i], matrix(l1_vector, nrow = 1), 
                                         matrix(l2_vector, nrow = 1), cut_samples[n_i,y], cut_samples[n_i,y+1])
          }
        }
        means <- apply(est_samples, 1, median)
        quantile_025 <- apply(est_samples, 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
        quantile_975 <- apply(est_samples, 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
        prediction[temp_index,] <- c(n_sample, y, round(c(means, quantile_025, quantile_975), digits = 4))
        temp_index <- temp_index + 1
      }
    }
    return(prediction)
  }
}
