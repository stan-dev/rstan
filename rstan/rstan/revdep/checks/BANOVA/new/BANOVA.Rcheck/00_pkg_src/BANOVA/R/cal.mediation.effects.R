# this function outputs effects of the mediator (sum of all coefficients related to 
# the mediator for each combination of moderators except the xvar)in the model sol_1
cal.mediation.effects <-
function (sol_1, est_matrix, n_sample, mediator, xvar = NA){
  table_list <- list()
  table_list_index <- 1
  model1_level1_var_matrix <- attr(attr(sol_1$mf1, 'terms'),'factors')
  model1_level1_var_dataClasses <- attr(attr(sol_1$mf1, 'terms'),'dataClasses')
  model1_level2_var_matrix <- attr(attr(sol_1$mf2, 'terms'),'factors')
  model1_level2_var_dataClasses <- attr(attr(sol_1$mf2, 'terms'),'dataClasses')
  
  # for each level of xvar, need to check xvar in level 1 or 2
  # find the mediator in the model and all other variables interacts with the mediator at the same level(moderators)
  # then find the coefficient of the mediator (for each combination of moderators) in model 1
  # algo:
  # e.g. (mediator, interactions) x (estimation matrix) x (all vars except xvar including other moderators), mediator(continuous, normal): value 1 since we only look at the coefficients
  # or  (all vars except xvar including other moderators) x (estimation matrix) x (mediator, interactions) if the mediator is at the btw. level
  # then find the coefficient of the xvar in model 2 (for each combination of moderators)
  # algo:
  # e.g. (xvar, interactions) x (estimation matrix) x (all vars except xvar including other moderators) 
  # or  (all vars except xvar including other moderators) x (estimation matrix) x (xvar, interactions) if the mediator is at the btw. level
  # Finally, join these two tables by common moderators
  # TODO: create a common function to deal with above two algos
  ### calculate mediation effects in model 1
  #if (attr(xvar, 'class') == 'numeric' || attr(xvar, 'class') == 'integer'){
  mediator_in_l1 <- mediator %in% rownames(model1_level1_var_matrix)
  mediator_in_l2 <- mediator %in% rownames(model1_level2_var_matrix)
  if (!mediator_in_l1 & !mediator_in_l2) stop(mediator," is not included in the model!")
  # mediator is at level 1
  if (mediator_in_l1){
    # find if there are moderators interacts with the mediator
    mediator_assign <- which (rownames(model1_level1_var_matrix) == mediator)
    interaction_list <- attr(sol_1$dMatrice$X, "interactions_num") # if the mediator is numeric 
    interaction_list_index <- attr(sol_1$dMatrice$X, "interactions_numeric_index")
    mediator_interaction_list <- list()
    mediator_interaction_list_index <- array()
    j = 1
    if (length(interaction_list) > 0){
      for (i in 1:length(interaction_list)){
        if (mediator_assign %in% interaction_list[[i]]){
          mediator_interaction_list[[j]] = interaction_list[[i]]
          mediator_interaction_list_index[j] = interaction_list_index[i]
          j = j + 1
        }
      }
    }
    interaction_list <- attr(sol_1$dMatrice$X, "interactions") # if the mediator is not numeric 
    interaction_list_index <- attr(sol_1$dMatrice$X, "interactions_index")
    if (length(interaction_list) > 0){
      for (i in 1:length(interaction_list)){
        if (mediator_assign %in% interaction_list[[i]]){
          mediator_interaction_list[[j]] = interaction_list[[i]]
          mediator_interaction_list_index[j] = interaction_list_index[i]
          j = j + 1
        }
      }
    }
    
    # for each interaction in mediator_interaction_list create a moderated mediation table (exclude main effects)
    # calculate l2 matrix first using all params but xvar, also using l2 formula
    # TOFIX: what about the case only intercept included
    l2_values <- attr(sol_1$dMatrice$Z, 'varValues')
    if (length(l2_values) == 0){
      # only intercept included
      l2_matrix <- model.matrix(~1)
      l2_matrix <- rbind(l2_matrix, c(1))
      attr(l2_matrix, "levels") <- l2_matrix
      if (sol_1$single_level){
        colnames(l2_matrix) <- c(" ")
      }
    }else{
      l2_matrix <- effect.matrix.mediator(interaction_factors = l2_values, matrix_formula=formula(attr(sol_1$mf2, 'terms')), xvar=xvar, contrast = sol_1$contrast)
    }
    if (length(mediator_interaction_list) > 0){
      l1_values <- attr(sol_1$dMatrice$X, 'varValues')
      mediator_interaction_effect_matrix <- list()
      for (i in 1:length(mediator_interaction_list)){
        # l1 matrix
        # TO check:
        # y var is also included in l1_values, interaction_list has considered this
        mediator_interaction_effect_matrix[[i]] <- effect.matrix.mediator(interaction_factors = l1_values[mediator_interaction_list[[i]]], mediator=mediator, xvar=xvar, contrast = sol_1$contrast)
        est_samples <- array(0, dim = c(nrow(mediator_interaction_effect_matrix[[i]]), nrow(l2_matrix), n_sample))
        for (n_s in 1:n_sample)
          est_samples[,,n_s] <- mediator_interaction_effect_matrix[[i]] %*% est_matrix[colnames(mediator_interaction_effect_matrix[[i]]), colnames(l2_matrix), n_s] %*% t(l2_matrix)
        table_m <- construct.table(est_samples, attr(mediator_interaction_effect_matrix[[i]], 'levels'), attr(l2_matrix, 'levels'))
        table_list[[table_list_index]] <- table_m
        table_list_index = table_list_index + 1
      }
      
    }else{
      l1_values <- attr(sol_1$dMatrice$X, 'varValues')
      # no interaction with the mediator in level 1, only select the mediator
      # TO check:
      # y var is also included in l1_values, interaction_list has considered this
      l1_matrix <- effect.matrix.mediator(l1_values[mediator_assign], mediator=mediator, contrast = sol_1$contrast)
      est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(l2_matrix), n_sample))
      for (n_s in 1:n_sample)
        est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(l2_matrix), n_s] %*% t(l2_matrix)
      
      table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'))
      table_list[[table_list_index]] <- table_m
      table_list_index = table_list_index + 1
    }
    
  }
  
  if (mediator_in_l2){
    # find if there are moderators interacts with the mediator
    mediator_assign <- which (rownames(model1_level2_var_matrix) == mediator)
    interaction_list <- attr(sol_1$dMatrice$Z, "interactions_num")
    interaction_list_index <- attr(sol_1$dMatrice$Z, "interactions_numeric_index")
    mediator_interaction_list <- list()
    mediator_interaction_list_index <- array()
    j = 1
    if (length(interaction_list) > 0){
      for (i in 1:length(interaction_list)){
        if (mediator_assign %in% interaction_list[[i]]){
          mediator_interaction_list[[j]] = interaction_list[[i]]
          mediator_interaction_list_index[j] = interaction_list_index[i]
          j = j + 1
        }
      }
    }
    interaction_list <- attr(sol_1$dMatrice$Z, "interactions") # if the mediator is not numeric 
    interaction_list_index <- attr(sol_1$dMatrice$Z, "interactions_index")
    if (length(interaction_list) > 0){
      for (i in 1:length(interaction_list)){
        if (mediator_assign %in% interaction_list[[i]]){
          mediator_interaction_list[[j]] = interaction_list[[i]]
          mediator_interaction_list_index[j] = interaction_list_index[i]
          j = j + 1
        }
      }
    }
    # for each interaction in mediator_interaction_list create a moderated mediation table (exclude main effects)
    # calculate l1 matrix first using all params but xvar, also using l1 formula
    l1_values <- attr(sol_1$dMatrice$X, 'varValues')
    # exclude y var
    l1_values <- l1_values[-1]
    if (length(l1_values) == 0){
      # only intercept included
      l1_matrix <- model.matrix(~1)
      l1_matrix <- rbind(l1_matrix, c(1))
      attr(l1_matrix, "levels") <- l1_matrix
    }else{
      l1_matrix <- effect.matrix.mediator(interaction_factors = l1_values, matrix_formula=formula(attr(sol_1$mf1, 'terms')), xvar=xvar, contrast = sol_1$contrast)
    }
    if (length(mediator_interaction_list) > 0){
      l2_values <- attr(sol_1$dMatrice$Z, 'varValues')
      mediator_interaction_effect_matrix <- list()
      for (i in 1:length(mediator_interaction_list)){
        # l2 matrix
        mediator_interaction_effect_matrix[[i]] <- effect.matrix.mediator(interaction_factors = l2_values[mediator_interaction_list[[i]]], mediator=mediator, xvar=xvar, contrast = sol_1$contrast)
        est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(mediator_interaction_effect_matrix[[i]]), n_sample))
        for (n_s in 1:n_sample){
          est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(mediator_interaction_effect_matrix[[i]]), n_s] %*% t(mediator_interaction_effect_matrix[[i]])
        }
        #TODO use a list to store
        table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(mediator_interaction_effect_matrix[[i]], 'levels'))
        table_list[[table_list_index]] <- table_m
        table_list_index = table_list_index + 1
      }
    }else{
      l2_values <- attr(sol_1$dMatrice$Z, 'varValues')
      # no interaction with the mediator in level 2, only select the mediator
      l2_matrix <- effect.matrix.mediator(l2_values[mediator_assign], mediator=mediator, contrast = sol_1$contrast)
      est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(l2_matrix), n_sample))
      for (n_s in 1:n_sample){
        est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(l2_matrix), n_s] %*% t(l2_matrix)
      }
      table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'))
      table_list[[table_list_index]] <- table_m
      table_list_index = table_list_index + 1
    }
  }
  return(table_list)
}

construct.table <- 
function (est_samples, row_name, col_name){
  means <- apply(est_samples, c(1,2), mean)
  quantile_025 <- apply(est_samples, c(1,2), quantile, probs = 0.025, type = 3, na.rm = FALSE)
  quantile_975 <- apply(est_samples, c(1,2), quantile, probs = 0.975, type = 3, na.rm = FALSE)
  table_m <- array(NA, dim = c(nrow(row_name) * nrow(col_name), ncol(row_name) + ncol(col_name) + 3), 
                 dimnames = list(rep("", nrow(row_name) * nrow(col_name)), c(colnames(row_name), colnames(col_name),'mean', '2.5%', '97.5%')))
  table_m_index <- array(NA, dim = c(nrow(row_name) * nrow(col_name), 2), 
                         dimnames = list(NULL, c('est_samples_row_index', 'est_samples_col_index')))
  for (k1 in 1:nrow(row_name)){
    temp <- ((k1-1) * nrow(col_name) + 1):((k1-1) * nrow(col_name) + nrow(col_name))
    table_m[temp, 1:ncol(row_name)] <- t(replicate(nrow(col_name), row_name[k1, ]))
    table_m[temp, (ncol(row_name) + 1) : (ncol(row_name) + ncol(col_name))] <- col_name
    table_m_index[temp,1] <- k1 
    table_m_index[temp,2] <- 1:nrow(col_name)
    table_m[temp, ncol(row_name) + ncol(col_name) + 1] <- round(means[k1,], digits = 4)
    table_m[temp, ncol(row_name) + ncol(col_name) + 2] <- pmin(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
    table_m[temp, ncol(row_name) + ncol(col_name) + 3] <- pmax(round(quantile_025[k1,], digits = 4), round(quantile_975[k1,], digits = 4))
  }
  # reorder table_m column names, sort values column by column, keep the order of table_m_index
  table_mediator <- list(table_m = table_m, index_name = table_m[, 1: (ncol(row_name) + ncol(col_name)), drop = F], index = table_m_index, samples = est_samples)
  return(table_mediator )
  
}
