est <-
function (f, est_matrix, n_sample, l1_matrix, l2_matrix, cutp0, cutp1){
  if (nrow(l1_matrix) == 1 & nrow(l2_matrix) == 1){
    est_samples_cutp0 <- array(0, dim = c(nrow(l2_matrix), n_sample))
    est_samples_cutp1 <- array(0, dim = c(nrow(l2_matrix), n_sample))
    est_samples <- array(0, dim = c(nrow(l2_matrix), n_sample))
    if (sum(cutp0 != cutp1) == 0){
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[, , n_s] %*% t(l2_matrix)
        est_samples[, n_s] <- 1 - f(temp)
      }
      
    }else if(sum(cutp1 != 0) == 0){
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[, , n_s] %*% t(l2_matrix) - cutp0[n_s]
        est_samples[, n_s] <- f(temp) 
      }
    }else{
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[, , n_s] %*% t(l2_matrix)
        est_samples_cutp0[, n_s] <- temp - cutp0[n_s]
        est_samples_cutp1[, n_s] <- temp - cutp1[n_s]
        est_samples[, n_s] <- f(est_samples_cutp0[, n_s]) - f(est_samples_cutp1[, n_s])
      }
    }
  }else if (nrow(l1_matrix) == 1 & nrow(l2_matrix) > 1){
    est_samples_cutp0 <- array(0, dim = c(nrow(l2_matrix), n_sample))
    est_samples_cutp1 <- array(0, dim = c(nrow(l2_matrix), n_sample))
    est_samples <- array(0, dim = c(nrow(l2_matrix), n_sample))
    if (sum(cutp0 != cutp1) == 0){
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[, colnames(l2_matrix), n_s] %*% t(l2_matrix)
        est_samples[, n_s] <- 1 - f(temp)
      }
      
    }else if(sum(cutp1 != 0) == 0){
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[, colnames(l2_matrix), n_s] %*% t(l2_matrix) - cutp0[n_s]
        est_samples[, n_s] <- f(temp) 
      }
    }else{
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[, colnames(l2_matrix), n_s] %*% t(l2_matrix)
        est_samples_cutp0[, n_s] <- temp - cutp0[n_s]
        est_samples_cutp1[, n_s] <- temp - cutp1[n_s]
        est_samples[, n_s] <- f(est_samples_cutp0[, n_s]) - f(est_samples_cutp1[, n_s])
      }
    }
  }else if (nrow(l2_matrix) == 1 & nrow(l1_matrix) > 1){
    est_samples_cutp0 <- array(0, dim = c(nrow(l1_matrix), n_sample))
    est_samples_cutp1 <- array(0, dim = c(nrow(l1_matrix), n_sample))
    est_samples <- array(0, dim = c(nrow(l1_matrix), n_sample))
    if (sum(cutp0 != cutp1) == 0){
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[colnames(l1_matrix), , n_s] %*% t(l2_matrix)
        est_samples[, n_s] <- 1 - f(temp)
      }
      
    }else if(sum(cutp1 != 0) == 0){
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[colnames(l1_matrix), , n_s] %*% t(l2_matrix) - cutp0[n_s]
        est_samples[, n_s] <- f(temp) 
      }
    }else{
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[colnames(l1_matrix), , n_s] %*% t(l2_matrix)
        est_samples_cutp0[, n_s] <- temp - cutp0[n_s]
        est_samples_cutp1[, n_s] <- temp - cutp1[n_s]
        est_samples[, n_s] <- f(est_samples_cutp0[, n_s]) - f(est_samples_cutp1[, n_s])
      }
    }
    
  }else{
    est_samples_cutp0 <- array(0, dim = c(nrow(l1_matrix), nrow(l2_matrix), n_sample))
    est_samples_cutp1 <- array(0, dim = c(nrow(l1_matrix), nrow(l2_matrix), n_sample))
    est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(l2_matrix), n_sample))
    if (sum(cutp0 != cutp1) == 0){
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(l2_matrix), n_s] %*% t(l2_matrix)
        est_samples[,, n_s] <- 1 - f(temp)
      }
      
    }else if(sum(cutp1 != 0) == 0){
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(l2_matrix), n_s] %*% t(l2_matrix) - cutp0[n_s]
        est_samples[,, n_s] <- f(temp) 
      }
    }else{
      for (n_s in 1:n_sample){
        temp <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(l2_matrix),n_s] %*% t(l2_matrix)
        est_samples_cutp0[,, n_s] <- temp - cutp0[n_s]
        est_samples_cutp1[,, n_s] <- temp - cutp1[n_s]
        est_samples[,, n_s] <- f(est_samples_cutp0[colnames(l1_matrix), colnames(l2_matrix), n_s]) - f(est_samples_cutp1[,, n_s])
      }
    }
  }
  return(est_samples)
}
