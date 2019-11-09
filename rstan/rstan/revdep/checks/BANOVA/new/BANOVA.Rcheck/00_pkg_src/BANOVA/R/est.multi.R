est.multi <-
function (est_matrix, n_sample, l1_matrix, l2_matrix){
  n_choice <- length(l1_matrix)
  exp_choice_linear <- list()
  choice_prob <- list()
  row_names <- dimnames(est_matrix)[[1]]
  for (n_c in 1: n_choice){
    if (nrow(l1_matrix[[n_c]]) == 1 & nrow(l2_matrix) == 1){
      est_samples <- array(0, dim = c(nrow(l2_matrix), n_sample))
      for (n_s in 1:n_sample){
        temp <- l1_matrix[[n_c]] %*% est_matrix[,,n_s] %*% t(l2_matrix)
        est_samples[, n_s] <- temp
      }
      exp_choice_linear[[n_c]] <- exp(est_samples)
    }else if (nrow(l1_matrix[[n_c]]) > 1 & nrow(l2_matrix) == 1){
      est_samples <- array(0, dim = c(nrow(l1_matrix[[n_c]]), n_sample))
      common_name <- intersect(colnames(l1_matrix[[n_c]]), row_names)
      l1_temp <- as.matrix(as.data.frame(l1_matrix[[n_c]])[common_name])
      for (n_s in 1:n_sample){
        # choice 1 intercept not in est_matrix
        temp <- l1_temp %*% est_matrix[common_name, , n_s] %*% t(l2_matrix)
        est_samples[, n_s] <- temp
      }
      exp_choice_linear[[n_c]] <- exp(est_samples)
    }else if (nrow(l1_matrix[[n_c]]) == 1 & nrow(l2_matrix) > 1){
      est_samples <- array(0, dim = c(nrow(l2_matrix), n_sample))
      for (n_s in 1:n_sample){
        temp <- l1_matrix[[n_c]] %*% est_matrix[, colnames(l2_matrix), n_s] %*% t(l2_matrix)
        est_samples[, n_s] <- temp
      }
      exp_choice_linear[[n_c]] <- exp(est_samples)
    }else{
      est_samples <- array(0, dim = c(nrow(l1_matrix[[n_c]]), nrow(l2_matrix), n_sample))
      for (n_s in 1:n_sample){
        common_name <- intersect(colnames(l1_matrix[[n_c]]), row_names)
        l1_temp <- as.matrix(as.data.frame(l1_matrix[[n_c]])[common_name])
        temp <- l1_temp %*% est_matrix[common_name, colnames(l2_matrix), n_s] %*% t(l2_matrix)
        est_samples[,, n_s] <- temp
      }
      exp_choice_linear[[n_c]] <- exp(est_samples)
    }
  }
  choice_prob_sum <- array(0 , dim = dim(exp_choice_linear[[1]]))
  for (n_c in 1:n_choice){
    choice_prob_sum <- choice_prob_sum + exp_choice_linear[[n_c]]
  }
  for (n_c in 1:n_choice){
    # print(exp_choice_linear[[n_c]])
    # print('sum:')
    # print(choice_prob_sum)
    # print('divid:')
    # print(exp_choice_linear[[n_c]]/choice_prob_sum)
    choice_prob[[n_c]] <-  exp_choice_linear[[n_c]]/choice_prob_sum
  }
  return(choice_prob)
}
