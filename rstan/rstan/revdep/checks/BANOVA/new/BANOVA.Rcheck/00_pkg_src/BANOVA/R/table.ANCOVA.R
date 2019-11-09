table.ANCOVA <-
  function(samples_l1_param, 
           X, 
           Z, 
           samples_l2_param, 
           y_val = NULL, 
           error = NULL, # for bernoulli binomial and multinomial case, sigle level, need to be fixed here
           multi = F, 
           n_cat = 0, 
           choice = 0,
           l1_error = 0 # variance of the level 1 error 
           ){
    eff_digits = 3
    if (length(attr(Z, 'varNames')) >= 1){
      if(!is.null(error)){ #for bernoulli binomial and multinomial case, single level
        num_l1_v <- ncol(X) # should be the same with num_l1
        num_l2_v <- length(attr(Z, 'varNames')) - 1# Not include intercept
        num_id <- nrow(Z)
        ancova_table <- data.frame(matrix(NA, nrow = num_l1_v, ncol = num_l2_v + 1 + 2)) # 1: intercept 
        rownames(ancova_table) <- colnames(X) #attr(X, 'varNames')
        colnames(ancova_table) <- c(attr(Z, 'varNames'), 'Residuals', 'Total') 
        effect_table <- data.frame(matrix(NA, nrow = num_l1_v, ncol = num_l2_v + 1)) # 1: intercept 
        rownames(effect_table) <- colnames(X) #attr(X, 'varNames')
        colnames(effect_table) <- c(attr(Z, 'varNames')) 
        X_assign <- attr(X, 'assign')
        Z_assign <- attr(Z, 'assign')
        if (multi)
          Z_factor_index <- 1:(num_l2_v+1)
        else
          Z_factor_index <- 0:num_l2_v # include intercept here
        n_f <- length(Z_factor_index)
        n_sample <- nrow(samples_l2_param)
        X_names <- colnames(X)
        Z_names <- colnames(Z)
        num_l1 <- length(X_assign)
        num_l2 <- length(Z_assign)
        est_matrix <- array(0 , dim = c(num_l1, num_l2, n_sample), dimnames = list(X_names, Z_names, NULL))
        for (i in 1:num_l1){
          for (j in 1:n_sample)
            est_matrix[i,,j] <- samples_l2_param[j,((i-1)*num_l2+1):((i-1)*num_l2+num_l2)]
        }
        for (i in 1:num_l1){
          # TODO: compute the error(finite-sample) for each factor
          factor_SS <- array(0, dim = c(1, n_f))
          factor_SS_numeric <- array(0, dim = c(1, n_f))
          effect_size <- array(0, dim = c(1, n_f))
          # get the new design matrix according to Z_assign and Z_factor_index
          newassign <- array(0, dim = 0) 
          index_of_factors <- array(0, dim = 0)
          for (i_z in 1:n_f){
            newassign <- c(newassign, Z_assign[Z_assign == Z_factor_index[i_z]])
            index_of_factors <- c(index_of_factors, which(Z_assign == Z_factor_index[i_z]))
          }
          newdesign_matrix <- Z[, index_of_factors]
          if (length(index_of_factors) > 1){
            if(!multi)
              if (qr(newdesign_matrix)$rank < length(index_of_factors)) stop ('Colinearity in design matrix, model is unidentified!')
          }
          pred_full <- array(0 , dim = c(num_id, n_sample))
          for (j in 1:n_sample){
            pred_full[ , j] <- est_matrix[i, , j] %*% t(newdesign_matrix) 
          }
          error <- array(error, dim = c(num_id, n_sample))
          var_error <- colSums(error) / num_id * (num_id - 1)
          e_error <- max(mean(var_error), 0)
          var_pred_full <- colSums(pred_full^2) / num_id * (num_id - 1)
          pred_full_center <- scale(pred_full, scale = F)
          var_pred_full_center <- colSums(pred_full_center^2) / num_id * (num_id - 1)
          #y_var <- colSums(y_center^2)#/(num_id - 1))
          #e_y_var <- mean(y_var)
          #y_square <- colSums(y^2)
          #ancova_table[i,'Residuals'] <- paste(format(round(e_error, digits = 4), nsmall = 4), " (", paste(pmax(round(quantile(var_error, c(0.05, 0.95)),digits = 2),0), collapse = ","), ")", sep = "")
          #ancova_table[i,'Total'] <- paste(format(round(e_y_var, digits = 4), nsmall = 4), " (", paste(round(quantile(y_var, c(0.05, 0.95)),digits = 2), collapse = ","), ")", sep = "")
          ancova_table[i,'Residuals'] <- format(round(e_error, digits = 4), nsmall = 4)
          if (multi){
            if (n_f >= (n_cat - 1) && choice > 0){
              # calculate the effect of the intercept
              factor_SS[choice] <- round(max(mean(var_pred_full - var_pred_full_center), 0), digits = 4)
              factor_SS_numeric[choice] <- factor_SS[choice]
              effect_size[choice] <- round(max(mean((var_pred_full - var_pred_full_center)/(var_pred_full - var_pred_full_center + var_error)), 0), digits = 4)
              #print(format(round(quantile(y_square - y_var, c(0.05, 0.95)),digits = 4),nsmall = 2))
              if (factor_SS[choice] > 0){
                #factor_SS[1] <- paste(format(factor_SS[1], nsmall = 4), " (", paste(pmax(round(quantile(y_square - y_var, c(0.05, 0.95)),digits = 2), 0), collapse = ","), ")", sep = "")
                factor_SS[choice] <- format(factor_SS[choice], nsmall = 4)
                #effect_size[1] <- paste(format(effect_size[1], nsmall = 4), " (", paste(pmax(round(quantile((y_square - y_var)/(y_square - y_var + var_error), c(0.05, 0.95)),digits = 2), 0), collapse = ","), ")", sep = "")
                effect_size[choice] <- paste(format(effect_size[choice], nsmall = 4), " (", paste(round(quantile((var_pred_full - var_pred_full_center)/(var_pred_full - var_pred_full_center + var_error), c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
              }else{
                #factor_SS[1] <- paste(format(factor_SS[1], nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
                factor_SS[choice] <- format(factor_SS[choice], nsmall = 4)
                #effect_size[1] <- paste(format(effect_size[1], nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
                effect_size[choice] <- paste(format(effect_size[choice], nsmall = 4), " (", paste(round(quantile((var_pred_full - var_pred_full_center)/(var_pred_full - var_pred_full_center + var_error), c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
              }
              #factor_SS[1] <- round(error, digits = 4)
              #paste(round(SS_full, digits = 5), ' (', round((SS_full)/SS_TO*100, digits = 2), '%)', sep="")
            }
            if (n_f >= n_cat){
              for (i_2 in n_cat:n_f){
                factor_index_i_2 <- which(newassign == Z_factor_index[i_2])
                design_matrix_i_2 <- newdesign_matrix
                design_matrix_i_2[, -factor_index_i_2] <- 0
                pred_i_2 <- array(0 , dim = c(num_id, n_sample))
                for (j in 1:n_sample){
                  pred_i_2[ , j] <- est_matrix[i, , j] %*% t(design_matrix_i_2) 
                }
                pred_i_2_center <- scale(pred_i_2, scale = F)
                error_pred_i_2_center <- colSums(pred_i_2_center^2) / num_id * (num_id - 1)
                var_error_i_2 <- error_pred_i_2_center + e_error
                var_error_i_2_type3 <- error_pred_i_2_center 
                e_error_i_2 <- max(mean(var_error_i_2_type3), 0)
                factor_SS_numeric[i_2] <- e_error_i_2
                effect_i_2 <- max(mean(var_error_i_2_type3/var_error_i_2), 0)
                if (e_error_i_2 > 0){
                  #factor_SS[i_2] <- paste(format(round(e_error_i_2, digits = 4), nsmall = 4), " (", paste(pmax(round(quantile(var_error_i_2_type3, c(0.05, 0.95)),digits = 2),0), collapse = ","), ")", sep = "")
                  factor_SS[i_2] <- format(round(e_error_i_2, digits = 4), nsmall = 4)
                  #effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(pmax(round(quantile(var_error_i_2_type3/var_error_i_2, c(0.05, 0.95)),digits = 2),0), collapse = ","), ")", sep = "")
                  effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(round(quantile(var_error_i_2_type3/var_error_i_2, c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
                }else{
                  #factor_SS[i_2] <- paste(format(round(e_error_i_2, digits = 4), nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
                  factor_SS[i_2] <- format(round(e_error_i_2, digits = 4), nsmall = 4)
                  #effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
                  effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(round(quantile(var_error_i_2_type3/var_error_i_2, c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
                }
                #paste(round(SS_full - SS_i, digits = 5), ' (',round((SS_full - SS_i)/SS_TO*100, digits = 2),'%)',sep = '')
              }
            }
            ancova_table[i, 1:n_f] = factor_SS
            ancova_table[i,'Total'] <- format(round(sum(factor_SS_numeric) + e_error, digits = 4), nsmall = 4)
            effect_table[i, 1:n_f] = effect_size
          }else{
            if (n_f >= 1){
              # calculate the effect of the intercept
              factor_SS[1] <- round(max(mean(var_pred_full - var_pred_full_center), 0), digits = 4)
              factor_SS_numeric[1] <- factor_SS[1]
              effect_size[1] <- round(max(mean((var_pred_full - var_pred_full_center)/(var_pred_full - var_pred_full_center + var_error)), 0), digits = 4)
              #print(format(round(quantile(y_square - y_var, c(0.05, 0.95)),digits = 4),nsmall = 2))
              if (factor_SS[1] > 0){
                #factor_SS[1] <- paste(format(factor_SS[1], nsmall = 4), " (", paste(pmax(round(quantile(y_square - y_var, c(0.05, 0.95)),digits = 2), 0), collapse = ","), ")", sep = "")
                factor_SS[1] <- format(factor_SS[1], nsmall = 4)
                #effect_size[1] <- paste(format(effect_size[1], nsmall = 4), " (", paste(pmax(round(quantile((y_square - y_var)/(y_square - y_var + var_error), c(0.05, 0.95)),digits = 2), 0), collapse = ","), ")", sep = "")
                effect_size[1] <- paste(format(effect_size[1], nsmall = 4), " (", paste(round(quantile((var_pred_full - var_pred_full_center)/(var_pred_full - var_pred_full_center + var_error), c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
              }else{
                #factor_SS[1] <- paste(format(factor_SS[1], nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
                factor_SS[1] <- format(factor_SS[1], nsmall = 4)
                #effect_size[1] <- paste(format(effect_size[1], nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
                effect_size[1] <- paste(format(effect_size[1], nsmall = 4), " (", paste(round(quantile((var_pred_full - var_pred_full_center)/(var_pred_full - var_pred_full_center + var_error), c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
              }
              #factor_SS[1] <- round(error, digits = 4)
              #paste(round(SS_full, digits = 5), ' (', round((SS_full)/SS_TO*100, digits = 2), '%)', sep="")
            }
            if (n_f >= 2){
              for (i_2 in 2:n_f){
                factor_index_i_2 <- which(newassign == Z_factor_index[i_2])
                design_matrix_i_2 <- newdesign_matrix
                design_matrix_i_2[, -factor_index_i_2] <- 0
                pred_i_2 <- array(0 , dim = c(num_id, n_sample))
                for (j in 1:n_sample){
                  pred_i_2[ , j] <- est_matrix[i, , j] %*% t(design_matrix_i_2) 
                }
                pred_i_2_center <- scale(pred_i_2, scale = F)
                error_pred_i_2_center <- colSums(pred_i_2_center^2) / num_id * (num_id - 1)
                var_error_i_2 <- error_pred_i_2_center + e_error
                var_error_i_2_type3 <- error_pred_i_2_center 
                e_error_i_2 <- max(mean(var_error_i_2_type3), 0)
                factor_SS_numeric[i_2] <- e_error_i_2
                effect_i_2 <- max(mean(var_error_i_2_type3/var_error_i_2), 0)
                if (e_error_i_2 > 0){
                  #factor_SS[i_2] <- paste(format(round(e_error_i_2, digits = 4), nsmall = 4), " (", paste(pmax(round(quantile(var_error_i_2_type3, c(0.05, 0.95)),digits = 2),0), collapse = ","), ")", sep = "")
                  factor_SS[i_2] <- format(round(e_error_i_2, digits = 4), nsmall = 4)
                  #effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(pmax(round(quantile(var_error_i_2_type3/var_error_i_2, c(0.05, 0.95)),digits = 2),0), collapse = ","), ")", sep = "")
                  effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(round(quantile(var_error_i_2_type3/var_error_i_2, c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
                }else{
                  #factor_SS[i_2] <- paste(format(round(e_error_i_2, digits = 4), nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
                  factor_SS[i_2] <- format(round(e_error_i_2, digits = 4), nsmall = 4)
                  #effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
                  effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(round(quantile(var_error_i_2_type3/var_error_i_2, c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
                }
                #paste(round(SS_full - SS_i, digits = 5), ' (',round((SS_full - SS_i)/SS_TO*100, digits = 2),'%)',sep = '')
              }
            }
            ancova_table[i, 1:n_f] = factor_SS
            ancova_table[i,'Total'] <- format(round(sum(factor_SS_numeric) + e_error, digits = 4), nsmall = 4)
            effect_table[i, 1:n_f] = effect_size
          }
        }
      }else{
        num_l1_v <- ncol(X) # should be the same with num_l1
        num_l2_v <- length(attr(Z, 'varNames')) - 1# Not include intercept
        num_id <- nrow(Z)
        ancova_table <- data.frame(matrix(NA, nrow = num_l1_v, ncol = num_l2_v + 1 + 2)) # 1: intercept 
        rownames(ancova_table) <- colnames(X) #attr(X, 'varNames')
        colnames(ancova_table) <- c(attr(Z, 'varNames'), 'Residuals', 'Total') 
        effect_table <- data.frame(matrix(NA, nrow = num_l1_v, ncol = num_l2_v + 1)) # 1: intercept 
        rownames(effect_table) <- colnames(X) #attr(X, 'varNames')
        colnames(effect_table) <- c(attr(Z, 'varNames')) 
        X_assign <- attr(X, 'assign')
        Z_assign <- attr(Z, 'assign')
        Z_factor_index <- 0:num_l2_v # include intercept here
        n_f <- length(Z_factor_index)
        n_sample <- nrow(samples_l2_param)
        X_names <- colnames(X)
        Z_names <- colnames(Z)
        
        num_l1 <- length(X_assign)
        num_l2 <- length(Z_assign)
        est_matrix <- array(0 , dim = c(num_l1, num_l2, n_sample), dimnames = list(X_names, Z_names, NULL))
        for (i in 1:num_l1){
          for (j in 1:n_sample)
            est_matrix[i,,j] <- samples_l2_param[j,((i-1)*num_l2+1):((i-1)*num_l2+num_l2)]
        }
        for (i in 1:num_l1){
          if (is.null(y_val))
            y <- t(as.matrix(samples_l1_param)[1:n_sample,((i-1)*num_id +1) : (i*num_id)])
          else
            y <- y_val
          # TODO: compute the error(finite-sample) for each factor
          factor_SS <- array(NA, dim = c(1, n_f))
          effect_size <- array(NA, dim = c(1, n_f))
          # get the new design matrix according to Z_assign and Z_factor_index
          newassign <- array(0, dim = 0) 
          index_of_factors <- array(0, dim = 0)
          for (i_z in 1:n_f){
            newassign <- c(newassign, Z_assign[Z_assign == Z_factor_index[i_z]])
            index_of_factors <- c(index_of_factors, which(Z_assign == Z_factor_index[i_z]))
          }
          newdesign_matrix <- Z[, index_of_factors]
          if (length(index_of_factors) > 1) 
            if (qr(newdesign_matrix)$rank < length(index_of_factors)) stop ('Colinearity in level 2 design matrix, model is unidentified!')
          
          pred_full <- array(0 , dim = c(num_id, n_sample))
          for (j in 1:n_sample){
            pred_full[ , j] <- est_matrix[i, , j] %*% t(newdesign_matrix) 
          }
          l1error <- array(l1_error, dim = c(num_id, n_sample))
          l1_var_error <- colSums(l1error) / num_id * (num_id - 1)
          if (!is.null(y_val)){
            error <- array(0 , dim = c(num_id, n_sample))
            for (p_col in 1:ncol(pred_full))
              error[, p_col] <- y_val - pred_full[, p_col]
          }else{
            error <- y - pred_full
          }
          error_center <- scale(error, scale = F)
          y_center <- scale(y, scale = F)
          var_error <- colSums(error_center^2) / num_id * (num_id - 1)
          e_error <- max(mean(var_error), 0)
          y_var <- colSums(y_center^2) / num_id * (num_id - 1)
          e_y_var <- mean(y_var)
          y_square <- colSums(y^2)
          #ancova_table[i,'Residuals'] <- paste(format(round(e_error, digits = 4), nsmall = 4), " (", paste(pmax(round(quantile(var_error, c(0.05, 0.95)),digits = 2),0), collapse = ","), ")", sep = "")
          #ancova_table[i,'Total'] <- paste(format(round(e_y_var, digits = 4), nsmall = 4), " (", paste(round(quantile(y_var, c(0.05, 0.95)),digits = 2), collapse = ","), ")", sep = "")
          ancova_table[i,'Residuals'] <- format(round(e_error, digits = 4), nsmall = 4)
          ancova_table[i,'Total'] <- format(round(mean(y_square), digits = 4), nsmall = 4)
          if (n_f >= 1){
            # calculate the effect of the intercept
            factor_SS[1] <- round(max(mean(y_square - y_var), 0), digits = 4)
            effect_size[1] <- round(max(mean((y_square - y_var)/(y_square - y_var + var_error + l1_var_error)), 0), digits = 4)
            #print(format(round(quantile(y_square - y_var, c(0.05, 0.95)),digits = 4),nsmall = 2))
            if (factor_SS[1] > 0){
              #factor_SS[1] <- paste(format(factor_SS[1], nsmall = 4), " (", paste(pmax(round(quantile(y_square - y_var, c(0.05, 0.95)),digits = 2), 0), collapse = ","), ")", sep = "")
              factor_SS[1] <- format(factor_SS[1], nsmall = 4)
              #effect_size[1] <- paste(format(effect_size[1], nsmall = 4), " (", paste(pmax(round(quantile((y_square - y_var)/(y_square - y_var + var_error), c(0.05, 0.95)),digits = 2), 0), collapse = ","), ")", sep = "")
              effect_size[1] <- paste(format(effect_size[1], nsmall = 4), " (", paste(round(quantile((y_square - y_var)/(y_square - y_var + var_error + l1_var_error), c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
            }else{
              #factor_SS[1] <- paste(format(factor_SS[1], nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
              factor_SS[1] <- format(factor_SS[1], nsmall = 4)
              #effect_size[1] <- paste(format(effect_size[1], nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
              effect_size[1] <- paste(format(effect_size[1], nsmall = 4), " (", paste(round(quantile((y_square - y_var)/(y_square - y_var + var_error + l1_var_error), c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
            }
            #factor_SS[1] <- round(error, digits = 4)
            #paste(round(SS_full, digits = 5), ' (', round((SS_full)/SS_TO*100, digits = 2), '%)', sep="")
          }
          if (n_f >= 2){
            for (i_2 in 2:n_f){
              factor_index_i_2 <- which(newassign == Z_factor_index[i_2])
              design_matrix_i_2 <- newdesign_matrix
              design_matrix_i_2[, -factor_index_i_2] <- 0
              pred_i_2 <- array(0 , dim = c(num_id, n_sample))
              for (j in 1:n_sample){
                pred_i_2[ , j] <- est_matrix[i, , j] %*% t(design_matrix_i_2) 
              }
              #print(head(design_matrix_i_2))
              error_i_2 <- error + pred_i_2
              error_i_2_center <- scale(error_i_2, scale = F)
              var_error_i_2 <- colSums(error_i_2_center^2) / num_id * (num_id - 1)
              var_error_i_2_type3 <- var_error_i_2 - var_error
              e_error_i_2 <- max(mean(var_error_i_2_type3), 0)
              effect_i_2 <- max(mean(var_error_i_2_type3/(var_error_i_2 + l1_var_error)), 0)
              if (e_error_i_2 > 0){
                #factor_SS[i_2] <- paste(format(round(e_error_i_2, digits = 4), nsmall = 4), " (", paste(pmax(round(quantile(var_error_i_2_type3, c(0.05, 0.95)),digits = 2),0), collapse = ","), ")", sep = "")
                factor_SS[i_2] <- format(round(e_error_i_2, digits = 4), nsmall = 4)
                #effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(pmax(round(quantile(var_error_i_2_type3/var_error_i_2, c(0.05, 0.95)),digits = 2),0), collapse = ","), ")", sep = "")
                effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(round(quantile(var_error_i_2_type3/(var_error_i_2 + l1_var_error), c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
              }else{
                #factor_SS[i_2] <- paste(format(round(e_error_i_2, digits = 4), nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
                factor_SS[i_2] <- format(round(e_error_i_2, digits = 4), nsmall = 4)
                #effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(c(0,0), collapse = ","), ")", sep = "")
                effect_size[i_2] <- paste(format(round(effect_i_2, digits = 4), nsmall = 4), " (", paste(round(quantile(var_error_i_2_type3/(var_error_i_2 + l1_var_error), c(0.025, 0.975)),digits = eff_digits), collapse = ","), ")", sep = "")
              }
              #paste(round(SS_full - SS_i, digits = 5), ' (',round((SS_full - SS_i)/SS_TO*100, digits = 2),'%)',sep = '')
            }
          }
          ancova_table[i, 1:n_f] = factor_SS
          effect_table[i, 1:n_f] = effect_size
        }
      }
      #return(data.frame(ancova_table))
      sol <- list(ancova_table = ancova_table, effect_table = effect_table)
      class(sol) <- 'ancova.effect'
      return(sol)
    }else
      return(NA)
  }