###
# In this version, the mediation analysis only includes one mediator
###
BANOVA.mediation <-
  function(sol_1, sol_2, xvar, mediator, individual = F){
    if(!(class(sol_1) %in% c('BANOVA', 'BANOVA.Normal', 'BANOVA.T', 'BANOVA.Poisson', 'BANOVA.Bern', 'BANOVA.Bin', 'BANOVA.ordMultinomial'))) stop('The Model is not supported yet')
    if(sol_1$model_name == 'BANOVA.Multinomial') stop('The Model is not supported yet')
    if(sol_2$model_name != 'BANOVA.Normal') stop('The mediator must follow the Normal distribution, use BANOVA Normal models instead.')
    
    X_names = colnames(sol_1$dMatrice$X)
    Z_names = colnames(sol_1$dMatrice$Z)
    X_assign = attr(sol_1$dMatrice$X, 'assign')
    Z_assign = attr(sol_1$dMatrice$Z, 'assign')
    num_l1 <- length(X_assign)
    num_l2 <- length(Z_assign)
    if (sol_1$single_level)
      samples_l2_param <- sol_1$samples_l1_param
    else
      samples_l2_param <- sol_1$samples_l2_param
    n_sample <- nrow(samples_l2_param)
    est_matrix <- array(0 , dim = c(num_l1, num_l2, n_sample), dimnames = list(X_names, Z_names, NULL))
    for (i in 1:num_l1){
      for (j in 1:n_sample)
        est_matrix[i,,j] <- samples_l2_param[j,((i-1)*num_l2+1):((i-1)*num_l2+num_l2)]
    }
    
    # mediator
    X_names_m = colnames(sol_2$dMatrice$X)
    Z_names_m = colnames(sol_2$dMatrice$Z)
    X_assign_m = attr(sol_2$dMatrice$X, 'assign')
    Z_assign_m = attr(sol_2$dMatrice$Z, 'assign')
    num_l1_m <- length(X_assign_m)
    num_l2_m <- length(Z_assign_m)
    if (sol_2$single_level)
      samples_l2_param_m <- sol_2$samples_l1_param
    else
      samples_l2_param_m <- sol_2$samples_l2_param
    n_sample_m <- nrow(samples_l2_param_m)
    est_matrix_m <- array(0 , dim = c(num_l1_m, num_l2_m, n_sample_m), dimnames = list(X_names_m, Z_names_m, NULL))
    for (i in 1:num_l1_m){
      for (j in 1:n_sample_m)
        est_matrix_m[i,,j] <- samples_l2_param_m[j,((i-1)*num_l2_m+1):((i-1)*num_l2_m+num_l2_m)]
    }
    sol <- list()
    
    if (individual){
      model1_level1_var_matrix <- attr(attr(sol_1$mf1, 'terms'),'factors')
      model1_level1_var_dataClasses <- attr(attr(sol_1$mf1, 'terms'),'dataClasses')
      model1_level2_var_matrix <- attr(attr(sol_1$mf2, 'terms'),'factors')
      model1_level2_var_dataClasses <- attr(attr(sol_1$mf2, 'terms'),'dataClasses')
      
      model2_level1_var_matrix <- attr(attr(sol_2$mf1, 'terms'),'factors')
      model2_level1_var_dataClasses <- attr(attr(sol_2$mf1, 'terms'),'dataClasses')
      model2_level2_var_matrix <- attr(attr(sol_2$mf2, 'terms'),'factors')
      model2_level2_var_dataClasses <- attr(attr(sol_2$mf2, 'terms'),'dataClasses')
      
      # used to extract the level 1 estimations fit_beta <- rstan::extract(out1$stan_fit, permuted = T)
      if (sol_1$single_level){
        stop("It seems to be a between-subject design, set individual = FALSE instead.")
      }
      fit_betas <- rstan::extract(sol_1$stan_fit, permuted = T)
      samples_l1_raw <- fit_betas$beta1
      #X_names = colnames(sol_1$dMatrice$X)
      samples_l1_individual <- aperm(samples_l1_raw, c(2,1,3)) # dimension: num_l1 x sample size x num_id
      dimnames(samples_l1_individual) <- list(X_names, NULL, NULL)
      
      # mediator
      fit_betas_m <- rstan::extract(sol_2$stan_fit, permuted = T)
      samples_l1_raw_m <- fit_betas_m$beta1
      #X_names_m = colnames(sol_2$dMatrice$X)
      samples_l1_individual_m <- aperm(samples_l1_raw_m, c(2,1,3)) # dimension: num_l1 x sample size x num_id
      dimnames(samples_l1_individual_m) <- list(X_names_m, NULL, NULL)
      
      # calc id map
      id_map = idmap(sol_1$old_id, sol_1$new_id)
      
      # calculate direct effect of xvar in model 1
      sol$dir_effects <- list()
      if (xvar %in% rownames(model1_level1_var_matrix)){
        sol$individual_direct <- list()
        direct_effects <- cal.mediation.effects.individual(sol_1, samples_l1_individual, xvar, mediator)
        for (i in 1:length(direct_effects)){
          # filter columns with only "1" or 1 (numeric), TODO: here 1 is hard coded, think about a better way to determine if it is numeric
          idx_to_rm <- c()
          for (j in 1:dim(direct_effects[[i]]$table_m[,,1,drop = F])[2]){
            if (all(direct_effects[[i]]$table_m[, j, 1] == '1') || all(direct_effects[[i]]$table_m[, j, 1] == 1))
              idx_to_rm <- c(idx_to_rm, j)
          }
          if(length(idx_to_rm) > 0)
            sol$dir_effects[[i]] <- direct_effects[[i]]$table_m[, -idx_to_rm, ,drop = F]
          else
            sol$dir_effects[[i]] <- direct_effects[[i]]$table_m
          sol$individual_direct[[i]] <- rabind(sol$dir_effects[[i]], id_map)
        }
      }else{
        direct_effects <- cal.mediation.effects(sol_1, est_matrix, n_sample, xvar, mediator)
        for (i in 1:length(direct_effects)){
          # filter columns with only "1" or 1 (numeric), TODO: here 1 is hard coded, think about a better way to determine if it is numeric
          idx_to_rm <- c()
          for (j in 1:ncol(direct_effects[[i]]$table_m)){
            if (all(direct_effects[[i]]$table_m[, j] == '1') || all(direct_effects[[i]]$table_m[, j] == 1))
              idx_to_rm <- c(idx_to_rm, j)
          }
          if(length(idx_to_rm) > 0)
            sol$dir_effects[[i]] <- direct_effects[[i]]$table_m[, -idx_to_rm]
          else
            sol$dir_effects[[i]] <- direct_effects[[i]]$table_m
          
        }
      }
      
      ## calculate effects of the mediator in model 1
      if (!(mediator %in% rownames(model1_level1_var_matrix))) stop("The mediator is between subjects, please set individual = FALSE.")
      mediator_l1_effects <- cal.mediation.effects.individual(sol_1, samples_l1_individual, mediator, xvar, is_mediator = T)
      sol$m1_effects <- list()
      for (i in 1:length(mediator_l1_effects)){
        # filter columns with only "1" or 1 (numeric), TODO: here 1 is hard coded, think about a better way to determine if it is numeric
        idx_to_rm <- c()
        for (j in 1:dim(mediator_l1_effects[[i]]$table_m[,,1, drop = F])[2]){
          if (all(mediator_l1_effects[[i]]$table_m[, j, 1] == '1') || all(mediator_l1_effects[[i]]$table_m[, j, 1] == 1))
            idx_to_rm <- c(idx_to_rm, j)
        }
        if(length(idx_to_rm) > 0)
          sol$m1_effects[[i]] <- mediator_l1_effects[[i]]$table_m[, -idx_to_rm, , drop = F]
        else
          sol$m1_effects[[i]] <- mediator_l1_effects[[i]]$table_m
        
      }
      
      # calculate effects of the xvar in model 2
      sol$m2_effects <- list()
      if (xvar %in% rownames(model2_level1_var_matrix)){ 
        mediator_xvar_effects <- cal.mediation.effects.individual(sol_2, samples_l1_individual_m, xvar)
        for (i in 1:length(mediator_xvar_effects)){
          mediator_xvar_effects[[i]]$table_m <- abind(array(1, dim = c(nrow(mediator_xvar_effects[[i]]$table_m), 1), dimnames = list(NULL, mediator)), mediator_xvar_effects[[i]]$table_m)
          mediator_xvar_effects[[i]]$index_name <- cbind(array(1, dim = c(nrow(mediator_xvar_effects[[i]]$index_name), 1), dimnames = list(NULL, mediator)), mediator_xvar_effects[[i]]$index_name)
          # filter columns with only "1" or 1 (numeric), TODO: here 1 is hard coded, think about a better way to determine if it is numeric
          idx_to_rm <- c()
          for (j in 1:dim(mediator_xvar_effects[[i]]$table_m[,,1, drop = F])[2]){
            if (all(mediator_xvar_effects[[i]]$table_m[, j, 1] == '1') || all(mediator_xvar_effects[[i]]$table_m[, j, 1] == 1))
              idx_to_rm <- c(idx_to_rm, j)
          }
          if(length(idx_to_rm) > 0)
            sol$m2_effects[[i]] <- mediator_xvar_effects[[i]]$table_m[, -idx_to_rm,, drop = F]
          else
            sol$m2_effects[[i]] <- mediator_xvar_effects[[i]]$table_m
          
        }
      }else if (xvar %in% rownames(model2_level2_var_matrix)){
        mediator_xvar_effects <- cal.mediation.effects(sol_2, est_matrix_m, n_sample_m, xvar)
        for (i in 1:length(mediator_xvar_effects)){
          mediator_xvar_effects[[i]]$table_m <- cbind(array(1, dim = c(nrow(mediator_xvar_effects[[i]]$table_m), 1), dimnames = list(NULL, mediator)), mediator_xvar_effects[[i]]$table_m)
          mediator_xvar_effects[[i]]$index_name <- cbind(array(1, dim = c(nrow(mediator_xvar_effects[[i]]$index_name), 1), dimnames = list(NULL, mediator)), mediator_xvar_effects[[i]]$index_name)
          # filter columns with only "1" or 1 (numeric), TODO: here 1 is hard coded, think about a better way to determine if it is numeric
          idx_to_rm <- c()
          for (j in 1:ncol(mediator_xvar_effects[[i]]$table_m)){
            if (all(mediator_xvar_effects[[i]]$table_m[, j] == '1') || all(mediator_xvar_effects[[i]]$table_m[, j] == 1))
              idx_to_rm <- c(idx_to_rm, j)
          }
          if(length(idx_to_rm) > 0)
            sol$m2_effects[[i]] <- mediator_xvar_effects[[i]]$table_m[, -idx_to_rm]
          else
            sol$m2_effects[[i]] <- mediator_xvar_effects[[i]]$table_m
          
        }
      }
      
      k <- 1
      sol$indir_effects <- list()
      sol$effect_size <- list()
      sol$individual_indirect <- list()
      for (i in 1:length(mediator_l1_effects))
        for (j in 1:length(mediator_xvar_effects)){
          comb_eff <- combine.effects.individual(mediator_l1_effects[[i]], mediator_xvar_effects[[j]], sol_1$tau_ySq, sol_1$data, mediator, id_map)
          indirect_effects <- comb_eff$table
          sol$effect_size[[k]] <- comb_eff$effect_size
          idx_to_rm <- c()
          for (j in 1:ncol(indirect_effects)){
            if (all(indirect_effects[, j, 1] == '1') || all(indirect_effects[, j, 1] == 1))
              idx_to_rm <- c(idx_to_rm, j)
          }
          if(length(idx_to_rm) > 0)
            sol$indir_effects[[k]] <- indirect_effects[, -idx_to_rm,, drop = F]
          else
            sol$indir_effects[[k]] <- indirect_effects
          sol$individual_indirect[[k]] <- rabind(sol$indir_effects[[k]])
          k <- k + 1
        }
      
      sol$xvar = xvar
      sol$mediator = mediator
      sol$individual = individual
      class(sol) <- 'BANOVA.mediation'
      return(sol)
      
    }else{
  
      #################
      if(0){ # temporarily replaced by the function cal.mediation.effects
      model1_level1_var_matrix <- attr(attr(sol_1$mf1, 'terms'),'factors')
      model1_level1_var_dataClasses <- attr(attr(sol_1$mf1, 'terms'),'dataClasses')
      model1_level2_var_matrix <- attr(attr(sol_1$mf2, 'terms'),'factors')
      model1_level2_var_dataClasses <- attr(attr(sol_1$mf2, 'terms'),'dataClasses')
      # find the direct effects
      # find corresponding names in X or Z for xvar, see floodlight analysis, then used in est_matrix
      xvar_in_l1 <- xvar %in% rownames(model1_level1_var_matrix)
      xvar_in_l2 <- xvar %in% rownames(model1_level2_var_matrix)
      if (!xvar_in_l1 & !xvar_in_l2) stop("xvar is not included in the model!")
      if (xvar_in_l1){
        attr(xvar, 'class') = model1_level1_var_dataClasses[xvar]
        xvar_index <- which(rownames(model1_level1_var_matrix) == xvar)
        xvar_index_assign <- which(model1_level1_var_matrix[xvar_index, ] == 1)
        xvar_related_names <- X_names[which(X_assign %in% xvar_index_assign)]
        for (xvar_name in xvar_related_names){
          print('Direct effects:')
          est_samples <- array(0, dim = c(n_sample))
          for (n_s in 1:n_sample){
            est_samples[n_s] <- est_matrix[xvar_name, 1, n_s]
          }
          direct_effec_mean <- mean(est_samples)
          quantiles <- quantile(est_samples, c(0.025, 0.975))
          tmp_output <- array(0, dim = c(1,3), dimnames = list(NULL, c('mean', '2.5%', '97.5%')))
          tmp_output[1,1] <- direct_effec_mean
          tmp_output[1,2:3] <- quantiles
          print(tmp_output)
        }
      }else{
        attr(xvar, 'class') = model1_level2_var_dataClasses[xvar]
        xvar_index <- which(rownames(model1_level2_var_matrix) == xvar)
        xvar_index_assign <- which(model1_level2_var_matrix[xvar_index, ] == 1)
        xvar_related_names <- Z_names[which(Z_assign %in% xvar_index_assign)]
        for (xvar_name in xvar_related_names){
          print('Direct effects:')
          est_samples <- array(0, dim = c(n_sample))
          for (n_s in 1:n_sample){
            est_samples[n_s] <- est_matrix[1, xvar_name, n_s]
          }
          direct_effec_mean <- mean(est_samples)
          quantiles <- quantile(est_samples, c(0.025, 0.975))
          tmp_output <- array(0, dim = c(1,3), dimnames = list(NULL, c('mean', '2.5%', '97.5%')))
          tmp_output[1,1] <- direct_effec_mean
          tmp_output[1,2:3] <- quantiles
          print(tmp_output)
        }
      }
      }
      ##################
      
      # calculate direct effect of xvar in model 1
      direct_effects <- cal.mediation.effects(sol_1, est_matrix, n_sample, xvar, mediator)
      sol$dir_effects <- list()
      for (i in 1:length(direct_effects)){
        # filter columns with only "1" or 1 (numeric), TODO: here 1 is hard coded, think about a better way to determine if it is numeric
        idx_to_rm <- c()
        for (j in 1:ncol(direct_effects[[i]]$table_m)){
          if (all(direct_effects[[i]]$table_m[, j] == '1') || all(direct_effects[[i]]$table_m[, j] == 1))
            idx_to_rm <- c(idx_to_rm, j)
        }
        if(length(idx_to_rm) > 0)
          sol$dir_effects[[i]] <- direct_effects[[i]]$table_m[, -idx_to_rm]
        else
          sol$dir_effects[[i]] <- direct_effects[[i]]$table_m
  
      }
      
      # calculate effects of the mediator in model 1
      mediator_l1_effects <- cal.mediation.effects(sol_1, est_matrix, n_sample, mediator, xvar)
      sol$m1_effects <- list()
      for (i in 1:length(mediator_l1_effects)){
        # filter columns with only "1" or 1 (numeric), TODO: here 1 is hard coded, think about a better way to determine if it is numeric
        idx_to_rm <- c()
        for (j in 1:ncol(mediator_l1_effects[[i]]$table_m)){
          if (all(mediator_l1_effects[[i]]$table_m[, j] == '1') || all(mediator_l1_effects[[i]]$table_m[, j] == 1))
            idx_to_rm <- c(idx_to_rm, j)
        }
        if(length(idx_to_rm) > 0)
          sol$m1_effects[[i]] <- mediator_l1_effects[[i]]$table_m[, -idx_to_rm]
        else
          sol$m1_effects[[i]] <- mediator_l1_effects[[i]]$table_m
        
      }
      
      # calculate effects of the xvar in model 2
      mediator_xvar_effects <- cal.mediation.effects(sol_2, est_matrix_m, n_sample_m, xvar)
      sol$m2_effects <- list()
      for (i in 1:length(mediator_xvar_effects)){
        mediator_xvar_effects[[i]]$table_m <- cbind(array(1, dim = c(nrow(mediator_xvar_effects[[i]]$table_m), 1), dimnames = list(NULL, mediator)), mediator_xvar_effects[[i]]$table_m)
        mediator_xvar_effects[[i]]$index_name <- cbind(array(1, dim = c(nrow(mediator_xvar_effects[[i]]$index_name), 1), dimnames = list(NULL, mediator)), mediator_xvar_effects[[i]]$index_name)
        # filter columns with only "1" or 1 (numeric), TODO: here 1 is hard coded, think about a better way to determine if it is numeric
        idx_to_rm <- c()
        for (j in 1:ncol(mediator_xvar_effects[[i]]$table_m)){
          if (all(mediator_xvar_effects[[i]]$table_m[, j] == '1') || all(mediator_xvar_effects[[i]]$table_m[, j] == 1))
            idx_to_rm <- c(idx_to_rm, j)
        }
        if(length(idx_to_rm) > 0)
          sol$m2_effects[[i]] <- mediator_xvar_effects[[i]]$table_m[, -idx_to_rm]
        else
          sol$m2_effects[[i]] <- mediator_xvar_effects[[i]]$table_m
  
      }
  
      k <- 1
      sol$indir_effects <- list()
      sol$effect_size <- list()
      for (i in 1:length(mediator_l1_effects))
        for (j in 1:length(mediator_xvar_effects)){
          comb_eff <- combine.effects (mediator_l1_effects[[i]], mediator_xvar_effects[[j]], sol_1$tau_ySq, sol_1$data, mediator)
          indirect_effects <- comb_eff$table
          sol$effect_size[[k]] <- comb_eff$effect_size
          idx_to_rm <- c()
          for (j in 1:ncol(indirect_effects)){
            if (all(indirect_effects[, j] == '1') || all(indirect_effects[, j] == 1))
              idx_to_rm <- c(idx_to_rm, j)
          }
          if(length(idx_to_rm) > 0)
              sol$indir_effects[[k]] <- indirect_effects[, -idx_to_rm]
          else
              sol$indir_effects[[k]] <- indirect_effects
          k <- k + 1
        }
      sol$xvar = xvar
      sol$mediator = mediator
      sol$individual = individual
      class(sol) <- 'BANOVA.mediation'
      return(sol)
    }
  }

combine.effects.individual <- function (mediator_l1_effects, mediator_xvar_effects, tau_ySq, data, mediator, id_map){
  num_id = dim(mediator_l1_effects$samples)[3]
  table_1_names <- mediator_l1_effects$index_name
  table_2_names <- mediator_xvar_effects$index_name
  # find common columns 
  temp_1 <- mediator_l1_effects$index
  colnames(temp_1) <- paste(colnames(temp_1), '.1', sep = "")
  table_1_names_index <- cbind(table_1_names, temp_1)
  temp_2 <- mediator_xvar_effects$index
  colnames(temp_2) <- paste(colnames(temp_2), '.2', sep = "")
  table_2_names_index <- cbind(table_2_names, temp_2)
  table_2_names_index.df <- table_2_names_index
  table_1_names_index.df <- table_1_names_index
  temp_table_index <- merge(table_2_names_index.df, table_1_names_index.df, by = intersect(colnames(table_1_names), colnames(table_2_names)), all.x = T)
  table_1_est_sample_index <- temp_table_index[,colnames(temp_1), drop = F]
  table_2_est_sample_index <- temp_table_index[,colnames(temp_2), drop = F]
  # standardize the names of this table, so that the output table looks consistant, direct vs indirect, e.g. sort the column names
  union_names <- union(colnames(table_1_names), colnames(table_2_names))
  union_names <- union_names[order(union_names)]
  result_table <- array('1', dim = c(nrow(temp_table_index), ncol(temp_table_index) - 2 + 3 + 1 + 2, num_id), dimnames = list(rep("",nrow(temp_table_index)), c(union_names, 'mean', '2.5%', '97.5%', 'p.value', 'id', 'effect size'), NULL))
  common_n_sample <- min(dim(mediator_l1_effects$samples)[2], dim(mediator_xvar_effects$samples)[2])
  result_table_sample <- array('1', dim = c(nrow(temp_table_index), ncol(temp_table_index) - 2 + common_n_sample, num_id), dimnames = list(rep("",nrow(temp_table_index)), c(union_names, paste('s_', 1:common_n_sample, sep = "")), NULL))
  effect_size <- rep("", num_id)
  for (i in 1:num_id){
    result_table[, 'id', i] <- id_map[i]
    for (nm in union_names){
      result_table[, nm, i] <- as.character(temp_table_index[[nm]])
      result_table_sample[, nm, i] <- as.character(temp_table_index[[nm]])
    }
    for (ind in 1:nrow(table_1_est_sample_index)){
      m_samples <- mediator_l1_effects$samples[table_1_est_sample_index[ind,1], 1:common_n_sample, i] * mediator_xvar_effects$samples[table_2_est_sample_index[ind,1], 1:common_n_sample, i]
      result_table[ind,'mean', i] <- round(mean(m_samples), 4)
      result_table[ind,c('2.5%', '97.5%'), i] <- round(quantile(m_samples, probs = c(0.025, 0.975)),4)
      result_table[ind,'p.value', i] <- ifelse(round(pValues(array(m_samples, dim = c(length(m_samples), 1))), 4) == 0, '<0.0001', round(pValues(array(m_samples, dim = c(length(m_samples), 1))), 4))
      result_table_sample[ind, paste('s_', 1:common_n_sample, sep = ""), i] <- m_samples
    }
    # compute effect size for the indirect effect
    if (mediator %in% union_names)
      union_names <- union_names[-which(union_names == mediator)]
    if ('(Intercept)' %in% union_names)
      union_names <- union_names[-which(union_names == '(Intercept)')]
    #remove numeric variable
    to_rm <- c()
    for (to_rm_ind in 1:length(union_names)){
      if (is.numeric(data[,union_names[to_rm_ind]])){
        to_rm <- c(to_rm, to_rm_ind)
      }
    }
    if (length(to_rm) > 0)
      union_names <- union_names[-to_rm]
    data_eff <- data[, union_names, drop = F]
    data_eff_sample <- merge(data_eff, result_table_sample[,,i], by = union_names, all.x = TRUE)
    data_eff_sample <- apply(data_eff_sample[, paste('s_', 1:common_n_sample, sep = "")], 2, as.character)
    data_eff_sample <- apply(data_eff_sample, 2, as.numeric)
    var_sample <- apply(data_eff_sample, 2, var)
    eff_sample <- var_sample/(var_sample + tau_ySq)
    effect_size[i] <- paste(round(mean(eff_sample), 3), " (", paste(round(quantile(eff_sample, probs = c(0.025, 0.975)),3), collapse = ','), ")", sep="")
    result_table[, 'effect size', i] <- effect_size[i]
    ##sort values column by column
    #result_table[,,i] <- data.frame(result_table[,,i], check.names=FALSE)
    #result_table <- result_table[do.call(order, result_table),,]
    
  }
  return(list(table = result_table, effect_size = effect_size))
}

combine.effects <- function (mediator_l1_effects, mediator_xvar_effects, tau_ySq, data, mediator){
  table_1_names <- mediator_l1_effects$index_name
  table_2_names <- mediator_xvar_effects$index_name
  # find common columns 
  temp_1 <- mediator_l1_effects$index
  colnames(temp_1) <- paste(colnames(temp_1), '.1', sep = "")
  table_1_names_index <- cbind(table_1_names, temp_1)
  temp_2 <- mediator_xvar_effects$index
  colnames(temp_2) <- paste(colnames(temp_2), '.2', sep = "")
  table_2_names_index <- cbind(table_2_names, temp_2)
  table_2_names_index.df <- table_2_names_index
  table_1_names_index.df <- table_1_names_index
  temp_table_index <- merge(table_2_names_index.df, table_1_names_index.df, by = intersect(colnames(table_1_names), colnames(table_2_names)), all.x = T)
  table_1_est_sample_index <- temp_table_index[,colnames(temp_1), drop = F]
  table_2_est_sample_index <- temp_table_index[,colnames(temp_2), drop = F]
  # standardize the names of this table, so that the output table looks consistant, direct vs indirect, e.g. sort the column names
  union_names <- union(colnames(table_1_names), colnames(table_2_names))
  union_names <- union_names[order(union_names)]
  result_table <- array('1', dim = c(nrow(temp_table_index), ncol(temp_table_index) - 4 + 3 + 1), dimnames = list(rep("",nrow(temp_table_index)), c(union_names, 'mean', '2.5%', '97.5%', 'p.value')))
  common_n_sample <- min(dim(mediator_l1_effects$samples)[3], dim(mediator_xvar_effects$samples)[3])
  result_table_sample <- array('1', dim = c(nrow(temp_table_index), ncol(temp_table_index) - 4 + common_n_sample), dimnames = list(rep("",nrow(temp_table_index)), c(union_names, paste('s_', 1:common_n_sample, sep = ""))))
  for (nm in union_names){
    result_table[, nm] <- as.character(temp_table_index[[nm]])
    result_table_sample[, nm] <- as.character(temp_table_index[[nm]])
  }
  for (ind in 1:nrow(table_1_est_sample_index)){
    m_samples <- mediator_l1_effects$samples[table_1_est_sample_index[ind,1], table_1_est_sample_index[ind,2], 1:common_n_sample] * mediator_xvar_effects$samples[table_2_est_sample_index[ind,1], table_2_est_sample_index[ind,2], 1:common_n_sample]
    result_table[ind,'mean'] <- round(mean(m_samples), 4)
    result_table[ind,c('2.5%', '97.5%')] <- round(quantile(m_samples, probs = c(0.025, 0.975)),4)
    result_table[ind,'p.value'] <- ifelse(round(pValues(array(m_samples, dim = c(length(m_samples), 1))), 4) == 0, '<0.0001', round(pValues(array(m_samples, dim = c(length(m_samples), 1))), 4))
    result_table_sample[ind, paste('s_', 1:common_n_sample, sep = "")] <- m_samples
  }
  # compute effect size for the indirect effect
  if (mediator %in% union_names)
    union_names <- union_names[-which(union_names == mediator)]
  if ('(Intercept)' %in% union_names)
    union_names <- union_names[-which(union_names == '(Intercept)')]
  # remove numeric variable
  to_rm <- c()
  for (i in 1:length(union_names)){
    if (is.numeric(data[,union_names[i]])){
      to_rm <- c(to_rm, i)
    }
  }
  if (length(to_rm) > 0)
    union_names <- union_names[-to_rm]
  data_eff <- data[, union_names, drop = F]
  data_eff_sample <- merge(data_eff, result_table_sample, by = union_names, all.x = TRUE)
  data_eff_sample <- apply(data_eff_sample[, paste('s_', 1:common_n_sample, sep = "")], 2, as.character)
  data_eff_sample <- apply(data_eff_sample, 2, as.numeric)
  var_sample <- apply(data_eff_sample, 2, var)
  eff_sample <- var_sample/(var_sample + tau_ySq)
  effect_size <- paste(round(mean(eff_sample), 3), " (", paste(round(quantile(eff_sample, probs = c(0.025, 0.975), na.rm = T),3), collapse = ','), ")", sep="")
  #sort values column by column
  result_table <- data.frame(result_table, check.names=FALSE)
  result_table <- result_table[do.call(order, result_table), ]
  return(list(table = result_table, effect_size = effect_size))
}

# a three dimensional array b cbind its first two dimension matrice with a matrix a (2 dim)
abind <- function(a, b){
  n_id <- dim(b)[3]
  n_r <- dim(b)[1]
  n_c <- dim(b)[2]
  res <- array(NA, dim = c(n_r, n_c+ncol(a), n_id), dimnames = list(NULL, c(colnames(a), dimnames(b)[[2]]), NULL))
  for (i in 1:n_id){
    res[,,i] <- cbind(a, b[,,i])
  }
  return(res)
}

# rbind all first two dimesion matrice for a 3 dimension matrix
rabind <- function(a, id_map = NULL){
  if(is.null(id_map)){
    res = array(0, dim = c(dim(a)[1], dim(a)[2],0))
    for (i in 1:dim(a)[3]){
      res <- rbind(res, a[,,i])
    }
  }else{
    # the last dimension is the new id
    res = array(0, dim = c(dim(a)[1], dim(a)[2] + 1,0))
    for (i in 1:dim(a)[3]){
      res <- rbind(res, cbind(a[,,i], id = id_map[i]))
    }
  }
  return(res)
}