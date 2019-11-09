multi.print.table.means <-
function (coeff_table, n_choice, samples_l2_param, X_names, X_assign = array(dim = 0), X_classes = character(0), Z_names, Z_assign = array(dim = 0), Z_classes = character(0), 
                               l1_values = list(), l1_interactions = list(), l1_interactions_index = array(dim = 0), 
                               l2_values = list(), l2_interactions = list(), l2_interactions_index = array(dim = 0), 
                               numeric_index_in_X, numeric_index_in_Z, single_level = F, contrast = NULL){
  
  sol_tables <- list()
  if (length(X_assign) == 1 && length(Z_assign) == 1) coeff_table <- matrix(coeff_table, nrow = 1) # the case that there is only one intercept 
  n_sample <- nrow(samples_l2_param)
  num_l1 <- length(X_assign)
  num_l2 <- length(Z_assign)
  est_matrix <- array(0 , dim = c(num_l1, num_l2, n_sample), dimnames = list(X_names, Z_names, NULL))
  for (i in 1:num_l1){
    for (j in 1:n_sample)
      est_matrix[i,,j] <- samples_l2_param[j,((i-1)*num_l2+1):((i-1)*num_l2+num_l2)]
  }
  
  if(single_level){
    cat('\n')
    cat('Table of probabilities for each category of the response\n')
    cat('--------------------------------------------------------\n')
    
    # Grand mean
    # cat('\n')
    # cat('Grand mean: \n')
    sol_tables[['Grand mean(each response): \n']] <- list()
    l1_v <- list()
    for (n_c in 1:n_choice){
      temp <- matrix(rep(0, num_l1), nrow = 1)
      if (n_c != 0)
        temp[n_c - 1] <- 1
      l1_v[[n_c]] <- temp
    }
    l2_v <- matrix(c(1, rep(0, num_l2 - 1)), nrow = 1)
    # gmean <- est.multi(est_matrix, n_sample, l1_v, l2_v)
    # for (n_c in 1:n_choice){
    #   cat('\n')
    #   cat('Choice:',n_c,'\n')
    #   cat(round(mean(gmean[[n_c]]), digits = 4), '\n')
    #   temp_title <- paste('Choice: ',n_c, sep ="")
    #   sol_tables[['Grand mean(each response): \n']][[temp_title]] <- as.table(matrix(round(quantile(gmean[[n_c]], probs = c(0.025, 0.975)), digits = 4), nrow = 1, ncol = 2, dimnames = list('',c('2.5%','97.5%'))))
    #   print(sol_tables[['Grand mean(each response): \n']][[temp_title]])
    # }
    
    # means of main effect in level 1 and 2 
    if (length(X_classes) != 0){
      l1_factors <- which(X_classes == 'factor')
      if (length(l1_factors) != 0){
        cat('\n')
        #cat('Means for factors at level 1: \n')
        sol_tables[['Means for factors at level 1: \n']] <- list()
        l1_matrix <- list()
        #l2_v <- matrix(c(1, rep(0, num_l2 - 1)), nrow = 1) # only look at the intercept in level 2
        for (i in 1:length(l1_factors)){
          l1_matrix[[i]] <- multi.effect.matrix.factor(n_choice, l1_values[[l1_factors[i]]], X_assign, l1_factors[i], numeric_index_in_X, contrast = contrast)
          # Compute median and quantile
          # print(dim(est_matrix))
          # print(est_matrix[,,1])
          # print(n_sample)
          # print(l2_v)
          est_samples <- est.multi(est_matrix, n_sample, l1_matrix[[i]], l2_v)
          cat('\nFactor:',attr(X_classes, 'names')[l1_factors[i]],'\n')
          for (n_c in 1: n_choice){
            cat('\n')
            cat('Choice:',n_c,'\n')
            means <- apply(est_samples[[n_c]], 1, mean)
            quantile_025 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
            quantile_975 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
            table <- matrix(NA, nrow = length(means), ncol = 4)
            table[, 2] <- means
            table[, 3] <- pmin(quantile_025, quantile_975)
            table[, 4] <- pmax(quantile_025, quantile_975)
            table <- round(table, digits = 4)
            table[, 1] <- attr(l1_matrix[[i]][[1]],'levels')
            colnames(table) <- c(attr(X_classes, 'names')[l1_factors[i]], 'mean', '2.5%', '97.5%')
            rownames(table) <- rep('', nrow(table))
            cat('\n')
            print(as.table(table))
            temp_title <- paste('Choice: ',n_c, sep ="")
            sol_tables[['Means for factors at level 1: \n']][[temp_title]] <- as.table(table)
          }
        }
      }
    }
    
    if (length(Z_classes) != 0){
      l2_factors <- which(Z_classes == 'factor')
      if (length(l2_factors) != 0){
        cat('\n')
        #cat('Means for factors at level 2: \n')
        sol_tables[['Means for factors at level 2: \n']] <- list()
        l2_matrix <- list()
        for (i in 1:length(l2_factors)){
          l2_matrix[[i]] <- multi.effect.matrix.factor(n_choice, l2_values[[l2_factors[i]]], Z_assign, l2_factors[i], numeric_index_in_Z, contrast = contrast)
          est_samples <- est.multi(est_matrix, n_sample, l1_v, l2_matrix[[i]])
          cat('\nFactor:',attr(Z_classes, 'names')[l2_factors[i]],'\n')
          for (n_c in 1: n_choice){
            cat('\n')
            cat('Choice:',n_c,'\n')
            means <- apply(est_samples[[n_c]], 1, mean)
            quantile_025 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
            quantile_975 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
            table <- matrix(NA, nrow = length(means), ncol = 4)
            table[, 2] <- means
            table[, 3] <- pmin(quantile_025, quantile_975)
            table[, 4] <- pmax(quantile_025, quantile_975)
            table <- round(table, digits = 4)
            table[, 1] <- attr(l2_matrix[[i]],'levels')
            colnames(table) <- c(attr(Z_classes, 'names')[l2_factors[i]], 'mean', '2.5%', '97.5%')
            rownames(table) <- rep('', nrow(table))
            cat('\n')
            print(as.table(table))
            temp_title <- paste('Choice: ',n_c, sep ="")
            sol_tables[['Means for factors at level 2: \n']][[temp_title]] <- as.table(table)
          }
        }
      }
    }
    
    # means of interactions between level 1 and 2
    if (length(X_classes) != 0 && length(Z_classes) != 0){
      if (length(l1_factors) != 0 && length(l2_factors) != 0){
        cat('\n')
        #cat('Means for interactions between level 1 and level 2 factors: \n')
        sol_tables[['Means for interactions between level 1 and level 2 factors: \n']] <- list()
        for (i in 1:length(l1_factors)){
          for (j in 1:length(l2_factors)){
            #means <- l1_matrix[[i]] %*% est_matrix_mean %*% t(l2_matrix[[j]])
            #quantile_025 <- l1_matrix[[i]] %*% est_matrix_025 %*% t(l2_matrix[[j]])
            #quantile_975 <- l1_matrix[[i]] %*% est_matrix_975 %*% t(l2_matrix[[j]])
            est_samples <- est.multi(est_matrix, n_sample, l1_matrix[[i]], l2_matrix[[j]])
            cat('\nFactor:',attr(X_classes, 'names')[l1_factors[i]],' ',attr(Z_classes, 'names')[l2_factors[j]],'\n')
            for (n_c in 1: n_choice){
              cat('\n')
              cat('Choice:',n_c,'\n')
              means <- apply(est_samples[[n_c]], c(1,2), mean)
              quantile_025 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.025, type = 3, na.rm = FALSE)
              quantile_975 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.975, type = 3, na.rm = FALSE)
              table <- matrix(NA, nrow = nrow(l1_matrix[[i]][[1]]) * nrow(l2_matrix[[j]]), ncol = 5)
              for (k1 in 1:nrow(l1_matrix[[i]][[1]])){
                temp <- ((k1-1) * nrow(l2_matrix[[j]]) + 1):((k1-1) * nrow(l2_matrix[[j]]) + nrow(l2_matrix[[j]]))
                table[temp, 3] <- round(means[k1,], digits = 4)
                table[temp, 4] <- pmin(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                table[temp, 5] <- pmax(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                table[temp, 1] <- rep(attr(l1_matrix[[i]][[1]],'levels')[k1], nrow(l2_matrix[[j]]))
                table[temp, 2] <- attr(l2_matrix[[j]],'levels')
              }
              colnames(table) <- c(attr(X_classes, 'names')[l1_factors[i]], attr(Z_classes, 'names')[l2_factors[j]], 'mean', '2.5%', '97.5%')
              rownames(table) <- rep('', nrow(table))
              cat('\n')
              print(as.table(table))
              temp_title <- paste('Choice: ',n_c, sep ="")
              sol_tables[['Means for interactions between level 1 and level 2 factors: \n']][[temp_title]] <- as.table(table)
            }
          }
        }          
      }    
    }
    
    # means of interactions in level 1 or 2
    if (length(l1_interactions) > 0){
      cat('\n')
      #cat('Means for interactions at level 1: \n')
      sol_tables[['Means for interactions at level 1: \n']] <- list()
      l1_inter_matrix <- list()
      index <- 1
      for (i in 1:length(l1_interactions)){
        temp1 <- l1_values[l1_interactions[[i]]]
        if (length(temp1) >= 2){
          #temp2 <- list()
          #for (j in 1:length(temp1))
          #temp2[[j]] <- levels(temp1[[j]])
          l1_inter_matrix[[index]] <- multi.effect.matrix.interaction(n_choice, interaction_factors = temp1, assign = X_assign, 
                                                                      l1_interactions[[i]], index_inter_factor = l1_interactions_index[i], 
                                                                      numeric_index_in_X, contrast = contrast) 
          est_samples <- est.multi(est_matrix, n_sample, l1_inter_matrix[[index]], l2_v)
          cat('\nFactor:',attr(X_classes, 'names')[l1_interactions[[i]] - 1],'\n')
          for (n_c in 1: n_choice){
            cat('\n')
            cat('Choice:',n_c,'\n')
            means <- apply(est_samples[[n_c]], 1, mean)
            quantile_025 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
            quantile_975 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
            table <- cbind(attr(l1_inter_matrix[[index]][[1]], 'levels'), 
                           round(means, digits = 4), 
                           pmin(round(quantile_025, digits = 4), round(quantile_975, digits = 4)), 
                           pmax(round(quantile_025, digits = 4), round(quantile_975, digits = 4)))
            colnames(table) <- c(attr(X_classes, 'names')[l1_interactions[[i]]], 'mean', '2.5%', '97.5%')
            rownames(table) <- rep('', nrow(table))
            cat('\n')
            print(as.table(table))
            temp_title <- paste('Choice: ', n_c, sep ="")
            sol_tables[['Means for interactions at level 1: \n']][[temp_title]] <- as.table(table)
          }
          
          # print the interaction between this interaction and level 2 factors if there exists
          if (length(Z_classes) != 0){
            l2_factors <- which(Z_classes == 'factor')
            if (length(l2_factors) != 0){
              cat('\n')
              #cat('Means for interactions between level 1 interactions and level 2 factors: \n')
              sol_tables[['Means for interactions between level 1 interactions and level 2 factors: \n']] <- list()
              for (j in 1:length(l2_factors)){
                #means <- l1_inter_matrix[[index]] %*% est_matrix_mean %*% t(l2_matrix[[j]])
                #quantile_025 <- l1_inter_matrix[[index]] %*% est_matrix_025 %*% t(l2_matrix[[j]])
                #quantile_975 <- l1_inter_matrix[[index]] %*% est_matrix_975 %*% t(l2_matrix[[j]])
                est_samples <- est.multi(est_matrix, n_sample, l1_inter_matrix[[index]], l2_matrix[[j]])
                cat('\nFactor:',attr(X_classes, 'names')[l1_interactions[[i]]],attr(Z_classes, 'names')[l2_factors[j]],'\n')
                for (n_c in 1: n_choice){
                  cat('\n')
                  cat('Choice:',n_c,'\n')
                  means <- apply(est_samples[[n_c]], c(1,2), mean)
                  quantile_025 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.025, type = 3, na.rm = FALSE)
                  quantile_975 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.975, type = 3, na.rm = FALSE)
                  
                  table <- matrix(NA, nrow = nrow(l1_inter_matrix[[i]][[1]]) * nrow(l2_matrix[[j]]), ncol = 6)
                  for (k1 in 1:nrow(l1_inter_matrix[[i]][[1]])){
                    temp <- ((k1-1) * nrow(l2_matrix[[j]]) + 1):((k1-1) * nrow(l2_matrix[[j]]) + nrow(l2_matrix[[j]]))
                    table[temp, 4] <- round(means[k1,], digits = 4)
                    table[temp, 5] <- pmin(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                    table[temp, 6] <- pmax(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                    table[temp, 1] <- attr(l1_inter_matrix[[i]][[1]],'levels')[k1,][1]
                    table[temp, 2] <- attr(l1_inter_matrix[[i]][[1]],'levels')[k1,][2]
                    table[temp, 3] <- attr(l2_matrix[[j]],'levels')
                  }
                  colnames(table) <- c(attr(X_classes, 'names')[l1_interactions[[i]] - 1], attr(Z_classes, 'names')[l2_factors[j]], 'mean', '2.5%', '97.5%')
                  rownames(table) <- rep('', nrow(table))
                  cat('\n')
                  print(as.table(table))
                  temp_title <- paste('Choice: ', n_c, sep ="")
                  sol_tables[['Means for interactions between level 1 interactions and level 2 factors: \n']][[temp_title]] <- as.table(table)
                }
              }
            }
          }
          index <- index + 1
        }
      }
    }
    
    if (length(l2_interactions) > 0){
      cat('\n')
      #cat('Means for interactions at level 2: \n')
      sol_tables[['Means for interactions at level 2: \n']] <- list()
      l2_inter_matrix <- list()
      index <- 1
      for (i in 1:length(l2_interactions)){
        temp1 <- l2_values[l2_interactions[[i]]]
        if (length(temp1) >= 2){
          #temp2 <- list()
          #for (j in 1:length(temp1))
          #temp2[[j]] <- levels(temp1[[j]])
          l2_inter_matrix[[index]] <- multi.effect.matrix.interaction(n_choice, interaction_factors = temp1, assign = Z_assign, 
                                                                l2_interactions[[i]], index_inter_factor = l2_interactions_index[i], 
                                                                numeric_index_in_Z, contrast = contrast) 
          est_samples <- est.multi(est_matrix, n_sample, l1_v, l2_inter_matrix[[index]])
          cat('\nFactor:',attr(Z_classes, 'names')[l2_interactions[[i]]],'\n')
          for (n_c in 1: n_choice){
            cat('\n')
            cat('Choice:',n_c,'\n')
            means <- apply(est_samples[[n_c]], 1, mean)
            quantile_025 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
            quantile_975 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
            table <- cbind(attr(l2_inter_matrix[[index]], 'levels'), 
                           matrix(round(means, digits = 4), ncol = 1), 
                           matrix(pmin(round(quantile_025, digits = 4), round(quantile_975, digits = 4)), ncol = 1), 
                           matrix(pmax(round(quantile_025, digits = 4), round(quantile_975, digits = 4)), ncol = 1))
            colnames(table) <- c(attr(Z_classes, 'names')[l2_interactions[[i]]], 'mean', '2.5%', '97.5%')
            rownames(table) <- rep('', nrow(table))
            cat('\n')
            print(as.table(table))
            temp_title <- paste('Choice: ', n_c, sep ="")
            sol_tables[['Means for interactions at level 2: \n']][[temp_title]] <- as.table(table)
          }
          
          # print the interaction between this interaction and level 1 factors if there exists
          if (length(X_classes) != 0){
            l1_factors <- which(X_classes == 'factor')
            if (length(l1_factors) != 0){
              cat('\n')
              cat('Means for interactions between level 2 interactions and level 1 factors: \n')
              sol_tables[['Means for interactions between level 2 interactions and level 1 factors: \n']] <- list()
              for (j in 1:length(l1_factors)){
                #means <- l1_matrix[[j]] %*% est_matrix_mean %*% t(l2_inter_matrix[[index]])
                #quantile_025 <- l1_matrix[[j]] %*% est_matrix_025 %*% t(l2_inter_matrix[[index]])
                #quantile_975 <- l1_matrix[[j]] %*% est_matrix_975 %*% t(l2_inter_matrix[[index]])
                est_samples <- est.multi(est_matrix, n_sample, l1_matrix[[j]], l2_inter_matrix[[index]])
                cat('\nFactor:',attr(X_classes, 'names')[l1_factors[[j]]],attr(Z_classes, 'names')[l2_interactions[[i]]], '\n')
                for (n_c in 1: n_choice){
                  cat('\n')
                  cat('Choice:',n_c,'\n')
                  means <- apply(est_samples[[n_c]], c(1,2), mean)
                  quantile_025 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.025, type = 3, na.rm = FALSE)
                  quantile_975 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.975, type = 3, na.rm = FALSE)
                  table <- matrix(NA, nrow = nrow(l1_matrix[[j]][[1]]) * nrow(l2_inter_matrix[[index]]), ncol = 6)
                  for (k1 in 1:nrow(l1_matrix[[j]][[1]])){
                    temp <- ((k1-1) * nrow(l2_inter_matrix[[index]]) + 1):((k1-1) * nrow(l2_inter_matrix[[index]]) + nrow(l2_inter_matrix[[index]]))
                    table[temp, 4] <- round(means[k1,], digits = 4)
                    table[temp, 5] <- pmin(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                    table[temp, 6] <- pmax(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                    table[temp, 1] <- attr(l1_matrix[[j]][[1]],'levels')[k1]
                    table[temp, 2:3] <- attr(l2_inter_matrix[[i]],'levels')
                  }
                  colnames(table) <- c(attr(X_classes, 'names')[l1_factors[[j]]], attr(Z_classes, 'names')[l2_interactions[[i]]], 'mean', '2.5%', '97.5%')
                  rownames(table) <- rep('', nrow(table))
                  cat('\n')
                  print(as.table(table))
                  temp_title <- paste('Choice: ', n_c, sep ="")
                  sol_tables[['Means for interactions between level 2 interactions and level 1 factors: \n']][[temp_title]] <- as.table(table)
                }
              }
            }
          }
          index <- index + 1
        }
      }
      
    }
  }else{
    # compute the overall table of means(prediction of y)
    cat('\n')
    cat('Table of means of the response\n')
    cat('------------------------------\n')
    # Grand mean
    cat('\n')
    cat('Grand mean: \n')
    l1_v <- list()
    for (n_c in 1:n_choice){
      temp <- matrix(rep(0, num_l1), nrow = 1)
      if (n_c != 0)
        temp[n_c - 1] <- 1
      l1_v[[n_c]] <- temp
    }
    l2_v <- matrix(c(1, rep(0, num_l2 - 1)), nrow = 1) # only look at the intercept in level 2
    gmean <- est.multi(est_matrix, n_sample, l1_v, l2_v)
    ggmean <- 0
    for (n_c in 1:n_choice){
      ggmean <- ggmean + n_c*gmean[[n_c]]
      
    }
    cat(round(mean(ggmean), digits = 4), '\n')
    sol_tables[['Grand mean: \n']] <- as.table(matrix(round(quantile(ggmean, probs = c(0.025, 0.975)), digits = 4), nrow = 1, ncol = 2, dimnames = list('',c('2.5%','97.5%'))))
    print(sol_tables[['Grand mean: \n']])
    
    # means by main effect in level 1 and 2 
    if (length(X_classes) != 0){
      l1_factors <- which(X_classes == 'factor')
      if (length(l1_factors) != 0){
        cat('\n')
        cat('Table of means for factors at level 1: \n')
        l1_matrix <- list()
        #l2_v <- matrix(c(1, rep(0, num_l2 - 1)), nrow = 1) # only look at the intercept in level 2
        for (i in 1:length(l1_factors)){
          l1_matrix[[i]] <- multi.effect.matrix.factor(n_choice, l1_values[[l1_factors[i]]], X_assign, l1_factors[i], numeric_index_in_X, contrast = contrast)
          # Compute median and quantile
          est_samples <- est.multi(est_matrix, n_sample, l1_matrix[[i]], l2_v)
          est_l1mean <- 0
          cat('\nFactor:',attr(X_classes, 'names')[l1_factors[i]],'\n')
          for (n_c in 1: n_choice){
            est_l1mean <- est_l1mean + n_c*est_samples[[n_c]]
          }
          means <- apply(est_l1mean, 1, mean)
          quantile_025 <- apply(est_l1mean, 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
          quantile_975 <- apply(est_l1mean, 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
          table <- matrix(NA, nrow = length(means), ncol = 4)
          table[, 2] <- means
          table[, 3] <- pmin(quantile_025, quantile_975)
          table[, 4] <- pmax(quantile_025, quantile_975)
          table <- round(table, digits = 4)
          table[, 1] <- attr(l1_matrix[[i]][[1]],'levels')
          colnames(table) <- c(attr(X_classes, 'names')[l1_factors[i]], 'mean', '2.5%', '97.5%')
          rownames(table) <- rep('', nrow(table))
          cat('\n')
          print(as.table(table))
          sol_tables[['Table of means for factors at level 1: \n']] <- as.table(table)
        }
      }
    }
    
    if (length(Z_classes) != 0){
      l2_factors <- which(Z_classes == 'factor')
      if (length(l2_factors) != 0){
        cat('\n')
        cat('Table of means of for factors at level 2: \n')
        l2_matrix <- list()
        for (i in 1:length(l2_factors)){
          l2_matrix[[i]] <- effect.matrix.factor(l2_values[[l2_factors[i]]], Z_assign, l2_factors[i], numeric_index_in_Z, contrast = contrast)
          est_samples <- est.multi(est_matrix, n_sample, l1_v, l2_matrix[[i]])
          est_l2mean <- 0
          cat('\nFactor:',attr(Z_classes, 'names')[l2_factors[i]],'\n')
          for (n_c in 1: n_choice){
            est_l2mean <- est_l2mean + n_c*est_samples[[n_c]]
          }
          means <- apply(est_l2mean, 1, mean)
          quantile_025 <- apply(est_l2mean, 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
          quantile_975 <- apply(est_l2mean, 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
          table <- matrix(NA, nrow = length(means), ncol = 4)
          table[, 2] <- means
          table[, 3] <- pmin(quantile_025, quantile_975)
          table[, 4] <- pmax(quantile_025, quantile_975)
          table <- round(table, digits = 4)
          table[, 1] <- attr(l2_matrix[[i]],'levels')
          colnames(table) <- c(attr(Z_classes, 'names')[l2_factors[i]], 'mean', '2.5%', '97.5%')
          rownames(table) <- rep('', nrow(table))
          cat('\n')
          print(as.table(table))
          sol_tables[['Table of means of for factors at level 2: \n']] <- as.table(table)
        }
      }
    }
    
    # means by interactions between level 1 and 2
    if (length(X_classes) != 0 && length(Z_classes) != 0){
      if (length(l1_factors) != 0 && length(l2_factors) != 0){
        cat('\n')
        cat('Table of means for interactions between level 1 and level 2 factors: \n')
        for (i in 1:length(l1_factors)){
          for (j in 1:length(l2_factors)){
            est_l1l2mean <- 0
            est_samples <- est.multi(est_matrix, n_sample, l1_matrix[[i]], l2_matrix[[j]])
            cat('\nFactor:',attr(X_classes, 'names')[l1_factors[i]],' ',attr(Z_classes, 'names')[l2_factors[j]],'\n')
            for (n_c in 1: n_choice){
              est_l1l2mean <- est_l1l2mean + n_c*est_samples[[n_c]]
            }
            
            means <- apply(est_l1l2mean, c(1,2), mean)
            quantile_025 <- apply(est_l1l2mean, c(1,2), quantile, probs = 0.025, type = 3, na.rm = FALSE)
            quantile_975 <- apply(est_l1l2mean, c(1,2), quantile, probs = 0.975, type = 3, na.rm = FALSE)
            table <- matrix(NA, nrow = nrow(l1_matrix[[i]][[1]]) * nrow(l2_matrix[[j]]), ncol = 5)
            for (k1 in 1:nrow(l1_matrix[[i]][[1]])){
              temp <- ((k1-1) * nrow(l2_matrix[[j]]) + 1):((k1-1) * nrow(l2_matrix[[j]]) + nrow(l2_matrix[[j]]))
              table[temp, 3] <- round(means[k1,], digits = 4)
              table[temp, 4] <- pmin(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
              table[temp, 5] <- pmax(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
              table[temp, 1] <- rep(attr(l1_matrix[[i]][[1]],'levels')[k1], nrow(l2_matrix[[j]]))
              table[temp, 2] <- attr(l2_matrix[[j]],'levels')
            }
            colnames(table) <- c(attr(X_classes, 'names')[l1_factors[i]], attr(Z_classes, 'names')[l2_factors[j]], 'mean', '2.5%', '97.5%')
            rownames(table) <- rep('', nrow(table))
            cat('\n')
            print(as.table(table)) 
            sol_tables[['Table of means for interactions between level 1 and level 2 factors: \n']] <- as.table(table)
          }
        }          
      }    
    }
    
    # means of interactions in level 1 or 2
    if (length(l1_interactions) > 0){
      cat('\n')
      cat('Table of means for interactions at level 1: \n')
      l1_inter_matrix <- list()
      index <- 1
      for (i in 1:length(l1_interactions)){
        temp1 <- l1_values[l1_interactions[[i]]]
        if (length(temp1) >= 2){
          #temp2 <- list()
          #for (j in 1:length(temp1))
            #temp2[[j]] <- levels(temp1[[j]])
          l1_inter_matrix[[index]] <- multi.effect.matrix.interaction(n_choice, interaction_factors = temp1, assign = X_assign, 
                                                                      l1_interactions[[i]], index_inter_factor = l1_interactions_index[i], 
                                                                      numeric_index_in_X, contrast = contrast)
          est_l1inmean <- 0
          est_samples <- est.multi(est_matrix, n_sample, l1_inter_matrix[[index]], l2_v)
          cat('\nFactor:',attr(X_classes, 'names')[l1_interactions[[i]] - 1],'\n')
          for (n_c in 1: n_choice){
            est_l1inmean <- est_l1inmean + n_c*est_samples[[n_c]]
          }
          
          means <- apply(est_l1inmean, 1, mean)
          quantile_025 <- apply(est_l1inmean, 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
          quantile_975 <- apply(est_l1inmean, 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
          table <- cbind(attr(l1_inter_matrix[[index]][[1]], 'levels'), 
                         round(means, digits = 4), 
                         pmin(round(quantile_025, digits = 4), round(quantile_975, digits = 4)), 
                         pmax(round(quantile_025, digits = 4), round(quantile_975, digits = 4)))
          colnames(table) <- c(attr(X_classes, 'names')[l1_interactions[[i]]], 'mean', '2.5%', '97.5%')
          rownames(table) <- rep('', nrow(table))
          cat('\n')
          print(as.table(table))
          sol_tables[['Table of means for interactions at level 1: \n']] <- as.table(table)
          
          # print the interaction between this interaction and level 2 factors if there exists
          if (length(Z_classes) != 0){
            l2_factors <- which(Z_classes == 'factor')
            if (length(l2_factors) != 0){
              cat('\n')
              cat('Table of means for interactions between level 1 interactions and level 2 factors: \n')
              for (j in 1:length(l2_factors)){
                est_l1inl2mean <- 0
                est_samples <- est.multi(est_matrix, n_sample, l1_inter_matrix[[index]], l2_matrix[[j]])
                cat('\nFactor:',attr(X_classes, 'names')[l1_interactions[[i]]],attr(Z_classes, 'names')[l2_factors[j]],'\n')
                for (n_c in 1: n_choice){
                  est_l1inl2mean <- est_l1inl2mean + n_c*est_samples[[n_c]]
                }
                
                means <- apply(est_l1inl2mean, c(1,2), mean)
                quantile_025 <- apply(est_l1inl2mean, c(1,2), quantile, probs = 0.025, type = 3, na.rm = FALSE)
                quantile_975 <- apply(est_l1inl2mean, c(1,2), quantile, probs = 0.975, type = 3, na.rm = FALSE)
                
                table <- matrix(NA, nrow = nrow(l1_inter_matrix[[i]][[1]]) * nrow(l2_matrix[[j]]), ncol = 6)
                for (k1 in 1:nrow(l1_inter_matrix[[i]][[1]])){
                  temp <- ((k1-1) * nrow(l2_matrix[[j]]) + 1):((k1-1) * nrow(l2_matrix[[j]]) + nrow(l2_matrix[[j]]))
                  table[temp, 4] <- round(means[k1,], digits = 4)
                  table[temp, 5] <- pmin(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                  table[temp, 6] <- pmax(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                  table[temp, 1] <- attr(l1_inter_matrix[[i]][[1]],'levels')[k1,][1]
                  table[temp, 2] <- attr(l1_inter_matrix[[i]][[1]],'levels')[k1,][2]
                  table[temp, 3] <- attr(l2_matrix[[j]],'levels')
                }
                colnames(table) <- c(attr(X_classes, 'names')[l1_interactions[[i]] - 1], attr(Z_classes, 'names')[l2_factors[j]], 'mean', '2.5%', '97.5%')
                rownames(table) <- rep('', nrow(table))
                cat('\n')
                print(as.table(table))
                sol_tables[['Table of means for interactions between level 1 interactions and level 2 factors: \n']] <- as.table(table)
              }
            }
          }
          index <- index + 1
        }
      }
    }
    
    if (length(l2_interactions) > 0){
      cat('\n')
      cat('Table of means for interactions at level 2: \n')
      l2_inter_matrix <- list()
      index <- 1
      for (i in 1:length(l2_interactions)){
        temp1 <- l2_values[l2_interactions[[i]]]
        if (length(temp1) >= 2){
          #temp2 <- list()
          #for (j in 1:length(temp1))
            #temp2[[j]] <- levels(temp1[[j]])
          l2_inter_matrix[[index]] <- effect.matrix.interaction(interaction_factors = temp1, assign = Z_assign, 
                                                                l2_interactions[[i]], index_inter_factor = l2_interactions_index[i], 
                                                                numeric_index_in_Z, contrast = contrast)
          est_l2inmean <- 0
          est_samples <- est.multi(est_matrix, n_sample, l1_v, l2_inter_matrix[[index]])
          cat('\nFactor:',attr(Z_classes, 'names')[l2_interactions[[i]]],'\n')
          for (n_c in 1: n_choice){
            est_l2inmean <- est_l2inmean + n_c*est_samples[[n_c]]
          }
          means <- apply(est_l2inmean, 1, mean)
          quantile_025 <- apply(est_l2inmean, 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
          quantile_975 <- apply(est_l2inmean, 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
          table <- cbind(attr(l2_inter_matrix[[index]], 'levels'), 
                         matrix(round(means, digits = 4), ncol = 1), 
                         matrix(pmin(round(quantile_025, digits = 4), round(quantile_975, digits = 4)), ncol = 1), 
                         matrix(pmax(round(quantile_025, digits = 4), round(quantile_975, digits = 4)), ncol = 1))
          colnames(table) <- c(attr(Z_classes, 'names')[l2_interactions[[i]]], 'mean', '2.5%', '97.5%')
          rownames(table) <- rep('', nrow(table))
          cat('\n')
          print(as.table(table))
          sol_tables[['Table of means for interactions at level 2: \n']] <- as.table(table)
          
          # print the interaction between this interaction and level 1 factors if there exists
          if (length(X_classes) != 0){
            l1_factors <- which(X_classes == 'factor')
            if (length(l1_factors) != 0){
              cat('\n')
              cat('Table of means for interactions between level 2 interactions and level 1 factors: \n')
              for (j in 1:length(l1_factors)){
                est_l2inl1mean <- 0
                est_samples <- est.multi(est_matrix, n_sample, l1_matrix[[j]], l2_inter_matrix[[index]])
                cat('\nFactor:',attr(X_classes, 'names')[l1_factors[[j]]],attr(Z_classes, 'names')[l2_interactions[[i]]], '\n')
                for (n_c in 1: n_choice){
                  est_l2inl1mean<- est_l2inl1mean  + n_c*est_samples[[n_c]]
                }
                means <- apply(est_l2inl1mean, c(1,2), mean)
                quantile_025 <- apply(est_l2inl1mean, c(1,2), quantile, probs = 0.025, type = 3, na.rm = FALSE)
                quantile_975 <- apply(est_l2inl1mean, c(1,2), quantile, probs = 0.975, type = 3, na.rm = FALSE)
                table <- matrix(NA, nrow = nrow(l1_matrix[[j]][[1]]) * nrow(l2_inter_matrix[[index]]), ncol = 6)
                for (k1 in 1:nrow(l1_matrix[[j]][[1]])){
                  temp <- ((k1-1) * nrow(l2_inter_matrix[[index]]) + 1):((k1-1) * nrow(l2_inter_matrix[[index]]) + nrow(l2_inter_matrix[[index]]))
                  table[temp, 4] <- round(means[k1,], digits = 4)
                  table[temp, 5] <- pmin(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                  table[temp, 6] <- pmax(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                  table[temp, 1] <- attr(l1_matrix[[j]][[1]],'levels')[k1]
                  table[temp, 2:3] <- attr(l2_inter_matrix[[i]],'levels')
                }
                colnames(table) <- c(attr(X_classes, 'names')[l1_factors[[j]]], attr(Z_classes, 'names')[l2_interactions[[i]]], 'mean', '2.5%', '97.5%')
                rownames(table) <- rep('', nrow(table))
                cat('\n')
                print(as.table(table))
                sol_tables[['Table of means for interactions between level 2 interactions and level 1 factors: \n']] <- as.table(table)
              }
            }
          }
          index <- index + 1
        }
      }
      
    }
    
    cat('\n')
    cat('Table of probabilities for each category of the response\n')
    cat('-------------------------------------------------------\n')
    
    # Grand mean
    cat('\n')
    cat('Grand mean: \n')
    sol_tables[['Grand mean(each response): \n']] <- list()
    l1_v <- list()
    for (n_c in 1:n_choice){
      temp <- matrix(rep(0, num_l1), nrow = 1)
      if (n_c != 0)
        temp[n_c - 1] <- 1
      l1_v[[n_c]] <- temp
    }
    l2_v <- matrix(c(1, rep(0, num_l2 - 1)), nrow = 1)
    gmean <- est.multi(est_matrix, n_sample, l1_v, l2_v)
    for (n_c in 1:n_choice){
      cat('\n')
      cat('Choice:',n_c,'\n')
      cat(round(mean(gmean[[n_c]]), digits = 4), '\n')
      temp_title <- paste('Choice: ',n_c, sep ="")
      sol_tables[['Grand mean(each response): \n']][[temp_title]] <- as.table(matrix(round(quantile(gmean[[n_c]], probs = c(0.025, 0.975)), digits = 4), nrow = 1, ncol = 2, dimnames = list('',c('2.5%','97.5%'))))
      print(sol_tables[['Grand mean(each response): \n']][[temp_title]])
    }
    
    # means of main effect in level 1 and 2 
    if (length(X_classes) != 0){
      l1_factors <- which(X_classes == 'factor')
      if (length(l1_factors) != 0){
        cat('\n')
        cat('Means for factors at level 1: \n')
        sol_tables[['Means for factors at level 1: \n']] <- list()
        l1_matrix <- list()
        #l2_v <- matrix(c(1, rep(0, num_l2 - 1)), nrow = 1) # only look at the intercept in level 2
        for (i in 1:length(l1_factors)){
          l1_matrix[[i]] <- multi.effect.matrix.factor(n_choice, l1_values[[l1_factors[i]]], X_assign, l1_factors[i], numeric_index_in_X, contrast = contrast)
          # Compute median and quantile
          est_samples <- est.multi(est_matrix, n_sample, l1_matrix[[i]], l2_v)
          cat('\nFactor:',attr(X_classes, 'names')[l1_factors[i]],'\n')
          for (n_c in 1: n_choice){
            cat('\n')
            cat('Choice:',n_c,'\n')
            means <- apply(est_samples[[n_c]], 1, mean)
            quantile_025 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
            quantile_975 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
            table <- matrix(NA, nrow = length(means), ncol = 4)
            table[, 2] <- means
            table[, 3] <- pmin(quantile_025, quantile_975)
            table[, 4] <- pmax(quantile_025, quantile_975)
            table <- round(table, digits = 4)
            table[, 1] <- attr(l1_matrix[[i]][[1]],'levels')
            colnames(table) <- c(attr(X_classes, 'names')[l1_factors[i]], 'mean', '2.5%', '97.5%')
            rownames(table) <- rep('', nrow(table))
            cat('\n')
            print(as.table(table))
            temp_title <- paste('Choice: ',n_c, sep ="")
            sol_tables[['Means for factors at level 1: \n']][[temp_title]] <- as.table(table)
          }
        }
      }
    }
    
    if (length(Z_classes) != 0){
      l2_factors <- which(Z_classes == 'factor')
      if (length(l2_factors) != 0){
        cat('\n')
        cat('Means for factors at level 2: \n')
        sol_tables[['Means for factors at level 2: \n']] <- list()
        l2_matrix <- list()
        for (i in 1:length(l2_factors)){
          l2_matrix[[i]] <- effect.matrix.factor(l2_values[[l2_factors[i]]], Z_assign, l2_factors[i], numeric_index_in_Z, contrast = contrast)
          est_samples <- est.multi(est_matrix, n_sample, l1_v, l2_matrix[[i]])
          cat('\nFactor:',attr(Z_classes, 'names')[l2_factors[i]],'\n')
          for (n_c in 1: n_choice){
            cat('\n')
            cat('Choice:',n_c,'\n')
            means <- apply(est_samples[[n_c]], 1, mean)
            quantile_025 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
            quantile_975 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
            table <- matrix(NA, nrow = length(means), ncol = 4)
            table[, 2] <- means
            table[, 3] <- pmin(quantile_025, quantile_975)
            table[, 4] <- pmax(quantile_025, quantile_975)
            table <- round(table, digits = 4)
            table[, 1] <- attr(l2_matrix[[i]],'levels')
            colnames(table) <- c(attr(Z_classes, 'names')[l2_factors[i]], 'mean', '2.5%', '97.5%')
            rownames(table) <- rep('', nrow(table))
            cat('\n')
            print(as.table(table))
            temp_title <- paste('Choice: ',n_c, sep ="")
            sol_tables[['Means for factors at level 2: \n']][[temp_title]] <- as.table(table)
          }
        }
      }
    }
    
    # means of interactions between level 1 and 2
    if (length(X_classes) != 0 && length(Z_classes) != 0){
      if (length(l1_factors) != 0 && length(l2_factors) != 0){
        cat('\n')
        cat('Means for interactions between level 1 and level 2 factors: \n')
        sol_tables[['Means for interactions between level 1 and level 2 factors: \n']] <- list()
        for (i in 1:length(l1_factors)){
          for (j in 1:length(l2_factors)){
            #means <- l1_matrix[[i]] %*% est_matrix_mean %*% t(l2_matrix[[j]])
            #quantile_025 <- l1_matrix[[i]] %*% est_matrix_025 %*% t(l2_matrix[[j]])
            #quantile_975 <- l1_matrix[[i]] %*% est_matrix_975 %*% t(l2_matrix[[j]])
            est_samples <- est.multi(est_matrix, n_sample, l1_matrix[[i]], l2_matrix[[j]])
            cat('\nFactor:',attr(X_classes, 'names')[l1_factors[i]],' ',attr(Z_classes, 'names')[l2_factors[j]],'\n')
            for (n_c in 1: n_choice){
              cat('\n')
              cat('Choice:',n_c,'\n')
              means <- apply(est_samples[[n_c]], c(1,2), mean)
              quantile_025 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.025, type = 3, na.rm = FALSE)
              quantile_975 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.975, type = 3, na.rm = FALSE)
              table <- matrix(NA, nrow = nrow(l1_matrix[[i]][[1]]) * nrow(l2_matrix[[j]]), ncol = 5)
              for (k1 in 1:nrow(l1_matrix[[i]][[1]])){
                temp <- ((k1-1) * nrow(l2_matrix[[j]]) + 1):((k1-1) * nrow(l2_matrix[[j]]) + nrow(l2_matrix[[j]]))
                table[temp, 3] <- round(means[k1,], digits = 4)
                table[temp, 4] <- pmin(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                table[temp, 5] <- pmax(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                table[temp, 1] <- rep(attr(l1_matrix[[i]][[1]],'levels')[k1], nrow(l2_matrix[[j]]))
                table[temp, 2] <- attr(l2_matrix[[j]],'levels')
              }
              colnames(table) <- c(attr(X_classes, 'names')[l1_factors[i]], attr(Z_classes, 'names')[l2_factors[j]], 'mean', '2.5%', '97.5%')
              rownames(table) <- rep('', nrow(table))
              cat('\n')
              print(as.table(table))
              temp_title <- paste('Choice: ',n_c, sep ="")
              sol_tables[['Means for interactions between level 1 and level 2 factors: \n']][[temp_title]] <- as.table(table)
            }
          }
        }          
      }    
    }
    
    # means of interactions in level 1 or 2
    if (length(l1_interactions) > 0){
      cat('\n')
      cat('Means for interactions at level 1: \n')
      sol_tables[['Means for interactions at level 1: \n']] <- list()
      l1_inter_matrix <- list()
      index <- 1
      for (i in 1:length(l1_interactions)){
        temp1 <- l1_values[l1_interactions[[i]]]
        if (length(temp1) >= 2){
          #temp2 <- list()
          #for (j in 1:length(temp1))
            #temp2[[j]] <- levels(temp1[[j]])
          l1_inter_matrix[[index]] <- multi.effect.matrix.interaction(n_choice, interaction_factors = temp1, assign = X_assign, 
                                                                l1_interactions[[i]], index_inter_factor = l1_interactions_index[i], 
                                                                numeric_index_in_X, contrast = contrast) 
          est_samples <- est.multi(est_matrix, n_sample, l1_inter_matrix[[index]], l2_v)
          cat('\nFactor:',attr(X_classes, 'names')[l1_interactions[[i]] - 1],'\n')
          for (n_c in 1: n_choice){
            cat('\n')
            cat('Choice:',n_c,'\n')
            means <- apply(est_samples[[n_c]], 1, mean)
            quantile_025 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
            quantile_975 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
            table <- cbind(attr(l1_inter_matrix[[index]][[1]], 'levels'), 
                           round(means, digits = 4), 
                           pmin(round(quantile_025, digits = 4), round(quantile_975, digits = 4)), 
                           pmax(round(quantile_025, digits = 4), round(quantile_975, digits = 4)))
            colnames(table) <- c(attr(X_classes, 'names')[l1_interactions[[i]]], 'mean', '2.5%', '97.5%')
            rownames(table) <- rep('', nrow(table))
            cat('\n')
            print(as.table(table))
            temp_title <- paste('Choice: ', n_c, sep ="")
            sol_tables[['Means for interactions at level 1: \n']][[temp_title]] <- as.table(table)
          }
          
          # print the interaction between this interaction and level 2 factors if there exists
          if (length(Z_classes) != 0){
            l2_factors <- which(Z_classes == 'factor')
            if (length(l2_factors) != 0){
              cat('\n')
              cat('Means for interactions between level 1 interactions and level 2 factors: \n')
              sol_tables[['Means for interactions between level 1 interactions and level 2 factors: \n']] <- list()
              for (j in 1:length(l2_factors)){
                #means <- l1_inter_matrix[[index]] %*% est_matrix_mean %*% t(l2_matrix[[j]])
                #quantile_025 <- l1_inter_matrix[[index]] %*% est_matrix_025 %*% t(l2_matrix[[j]])
                #quantile_975 <- l1_inter_matrix[[index]] %*% est_matrix_975 %*% t(l2_matrix[[j]])
                est_samples <- est.multi(est_matrix, n_sample, l1_inter_matrix[[index]], l2_matrix[[j]])
                cat('\nFactor:',attr(X_classes, 'names')[l1_interactions[[i]]],attr(Z_classes, 'names')[l2_factors[j]],'\n')
                for (n_c in 1: n_choice){
                  cat('\n')
                  cat('Choice:',n_c,'\n')
                  means <- apply(est_samples[[n_c]], c(1,2), mean)
                  quantile_025 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.025, type = 3, na.rm = FALSE)
                  quantile_975 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.975, type = 3, na.rm = FALSE)
                  
                  table <- matrix(NA, nrow = nrow(l1_inter_matrix[[i]][[1]]) * nrow(l2_matrix[[j]]), ncol = 6)
                  for (k1 in 1:nrow(l1_inter_matrix[[i]][[1]])){
                    temp <- ((k1-1) * nrow(l2_matrix[[j]]) + 1):((k1-1) * nrow(l2_matrix[[j]]) + nrow(l2_matrix[[j]]))
                    table[temp, 4] <- round(means[k1,], digits = 4)
                    table[temp, 5] <- pmin(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                    table[temp, 6] <- pmax(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                    table[temp, 1] <- attr(l1_inter_matrix[[i]][[1]],'levels')[k1,][1]
                    table[temp, 2] <- attr(l1_inter_matrix[[i]][[1]],'levels')[k1,][2]
                    table[temp, 3] <- attr(l2_matrix[[j]],'levels')
                  }
                  colnames(table) <- c(attr(X_classes, 'names')[l1_interactions[[i]] - 1], attr(Z_classes, 'names')[l2_factors[j]], 'mean', '2.5%', '97.5%')
                  rownames(table) <- rep('', nrow(table))
                  cat('\n')
                  print(as.table(table))
                  temp_title <- paste('Choice: ', n_c, sep ="")
                  sol_tables[['Means for interactions between level 1 interactions and level 2 factors: \n']][[temp_title]] <- as.table(table)
                }
              }
            }
          }
          index <- index + 1
        }
      }
    }
    
    if (length(l2_interactions) > 0){
      cat('\n')
      cat('Means for interactions at level 2: \n')
      sol_tables[['Means for interactions at level 2: \n']] <- list()
      l2_inter_matrix <- list()
      index <- 1
      for (i in 1:length(l2_interactions)){
        temp1 <- l2_values[l2_interactions[[i]]]
        if (length(temp1) >= 2){
          #temp2 <- list()
          #for (j in 1:length(temp1))
            #temp2[[j]] <- levels(temp1[[j]])
          l2_inter_matrix[[index]] <- effect.matrix.interaction(interaction_factors = temp1, assign = Z_assign, 
                                                                l2_interactions[[i]], index_inter_factor = l2_interactions_index[i], 
                                                                numeric_index_in_Z, contrast = contrast) 
          est_samples <- est.multi(est_matrix, n_sample, l1_v, l2_inter_matrix[[index]])
          cat('\nFactor:',attr(Z_classes, 'names')[l2_interactions[[i]]],'\n')
          for (n_c in 1: n_choice){
            cat('\n')
            cat('Choice:',n_c,'\n')
            means <- apply(est_samples[[n_c]], 1, mean)
            quantile_025 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.025, type = 3, na.rm = FALSE)
            quantile_975 <- apply(est_samples[[n_c]], 1, quantile, probs = 0.975, type = 3, na.rm = FALSE)
            table <- cbind(attr(l2_inter_matrix[[index]], 'levels'), 
                           matrix(round(means, digits = 4), ncol = 1), 
                           matrix(pmin(round(quantile_025, digits = 4), round(quantile_975, digits = 4)), ncol = 1), 
                           matrix(pmax(round(quantile_025, digits = 4), round(quantile_975, digits = 4)), ncol = 1))
            colnames(table) <- c(attr(Z_classes, 'names')[l2_interactions[[i]]], 'mean', '2.5%', '97.5%')
            rownames(table) <- rep('', nrow(table))
            cat('\n')
            print(as.table(table))
            temp_title <- paste('Choice: ', n_c, sep ="")
            sol_tables[['Means for interactions at level 2: \n']][[temp_title]] <- as.table(table)
          }
          
          # print the interaction between this interaction and level 1 factors if there exists
          if (length(X_classes) != 0){
            l1_factors <- which(X_classes == 'factor')
            if (length(l1_factors) != 0){
              cat('\n')
              cat('Means for interactions between level 2 interactions and level 1 factors: \n')
              sol_tables[['Means for interactions between level 2 interactions and level 1 factors: \n']] <- list()
              for (j in 1:length(l1_factors)){
                #means <- l1_matrix[[j]] %*% est_matrix_mean %*% t(l2_inter_matrix[[index]])
                #quantile_025 <- l1_matrix[[j]] %*% est_matrix_025 %*% t(l2_inter_matrix[[index]])
                #quantile_975 <- l1_matrix[[j]] %*% est_matrix_975 %*% t(l2_inter_matrix[[index]])
                est_samples <- est.multi(est_matrix, n_sample, l1_matrix[[j]], l2_inter_matrix[[index]])
                cat('\nFactor:',attr(X_classes, 'names')[l1_factors[[j]]],attr(Z_classes, 'names')[l2_interactions[[i]]], '\n')
                for (n_c in 1: n_choice){
                  cat('\n')
                  cat('Choice:',n_c,'\n')
                  means <- apply(est_samples[[n_c]], c(1,2), mean)
                  quantile_025 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.025, type = 3, na.rm = FALSE)
                  quantile_975 <- apply(est_samples[[n_c]], c(1,2), quantile, probs = 0.975, type = 3, na.rm = FALSE)
                  table <- matrix(NA, nrow = nrow(l1_matrix[[j]][[1]]) * nrow(l2_inter_matrix[[index]]), ncol = 6)
                  for (k1 in 1:nrow(l1_matrix[[j]][[1]])){
                    temp <- ((k1-1) * nrow(l2_inter_matrix[[index]]) + 1):((k1-1) * nrow(l2_inter_matrix[[index]]) + nrow(l2_inter_matrix[[index]]))
                    table[temp, 4] <- round(means[k1,], digits = 4)
                    table[temp, 5] <- pmin(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                    table[temp, 6] <- pmax(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
                    table[temp, 1] <- attr(l1_matrix[[j]][[1]],'levels')[k1]
                    table[temp, 2:3] <- attr(l2_inter_matrix[[i]],'levels')
                  }
                  colnames(table) <- c(attr(X_classes, 'names')[l1_factors[[j]]], attr(Z_classes, 'names')[l2_interactions[[i]]], 'mean', '2.5%', '97.5%')
                  rownames(table) <- rep('', nrow(table))
                  cat('\n')
                  print(as.table(table))
                  temp_title <- paste('Choice: ', n_c, sep ="")
                  sol_tables[['Means for interactions between level 2 interactions and level 1 factors: \n']][[temp_title]] <- as.table(table)
                }
              }
            }
          }
          index <- index + 1
        }
      }
      
    }
  }
  return(sol_tables)
}
