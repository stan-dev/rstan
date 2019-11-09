# this function outputs effects of the mediator (sum of all coefficients related to 
# the mediator for each combination of moderators except the xvar)in the model sol_1 at the individual level
cal.mediation.effects.individual <- function (sol_1, 
                                              est_matrix, 
                                              mediator, 
                                              xvar = NA,
                                              is_mediator = F){
  table_list <- list()
  table_list_index <- 1
  model1_level1_var_matrix <- attr(attr(sol_1$mf1, 'terms'),'factors')
  model1_level1_var_dataClasses <- attr(attr(sol_1$mf1, 'terms'),'dataClasses')
  model1_level2_var_matrix <- attr(attr(sol_1$mf2, 'terms'),'factors')
  model1_level2_var_dataClasses <- attr(attr(sol_1$mf2, 'terms'),'dataClasses')
  num_l1 = dim(est_matrix)[1]
  n_sample = dim(est_matrix)[2]
  num_id = dim(est_matrix)[3]
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
    
    if (length(mediator_interaction_list) > 0){
      l1_values <- attr(sol_1$dMatrice$X, 'varValues')
      mediator_interaction_effect_matrix <- list()
      for (i in 1:length(mediator_interaction_list)){
        # l1 matrix
        # TO check:
        # y var is also included in l1_values, interaction_list has considered this
        mediator_interaction_effect_matrix[[i]] <- effect.matrix.mediator(interaction_factors = l1_values[mediator_interaction_list[[i]]], mediator=mediator, xvar=xvar, contrast = sol_1$contrast)
        est_samples <- array(0, dim = c(nrow(mediator_interaction_effect_matrix[[i]]), n_sample, num_id))
        for (n_id in 1:num_id)
          est_samples[,,n_id] <- mediator_interaction_effect_matrix[[i]] %*% est_matrix[colnames(mediator_interaction_effect_matrix[[i]]),, n_id]
        table_m <- construct.table.individual(est_samples, attr(mediator_interaction_effect_matrix[[i]], 'levels'))
        table_list[[table_list_index]] <- table_m
        table_list_index = table_list_index + 1
      }
      
    }else{
      l1_values <- attr(sol_1$dMatrice$X, 'varValues')
      # no interaction with the mediator in level 1, only select the mediator
      # TO check:
      # y var is also included in l1_values, interaction_list has considered this
      l1_matrix <- effect.matrix.mediator(l1_values[mediator_assign], mediator=mediator, contrast = sol_1$contrast)
      est_samples <- array(0, dim = c(nrow(l1_matrix), n_sample, num_id))
      for (n_id in 1:num_id)
        est_samples[,,n_id] <- l1_matrix %*% est_matrix[colnames(l1_matrix),, n_id]
      
      table_m <- construct.table.individual(est_samples, attr(l1_matrix, 'levels'))
      table_list[[table_list_index]] <- table_m
      table_list_index = table_list_index + 1
    }
    
  }
  
  if (mediator_in_l2){ # xvar in level 2
    if (is_mediator){
      stop("The mediator is between subjects, please set individual = FALSE.")
    }

  }
  return(table_list)
}

construct.table.individual <- 
  function (est_samples, row_name){
    num_id = dim(est_samples)[3]
    means <- apply(est_samples, c(1,3), mean)
    quantile_025 <- apply(est_samples, c(1,3), quantile, probs = 0.025, type = 3, na.rm = FALSE)
    quantile_975 <- apply(est_samples, c(1,3), quantile, probs = 0.975, type = 3, na.rm = FALSE)
    table_m <- array(NA, dim = c(nrow(row_name), ncol(row_name) + 3, num_id), 
                     dimnames = list(rep("", nrow(row_name)), c(colnames(row_name), 'mean', '2.5%', '97.5%'), NULL))
    table_m_index <- array(NA, dim = c(nrow(row_name), 1), 
                           dimnames = list(NULL, c('est_samples_row_index')))
    for (n_id in 1:num_id){
      for (k1 in 1:nrow(row_name)){
        temp <- k1
        table_m[temp, 1:ncol(row_name), n_id] <- row_name[k1, ]
        #table_m[temp, (ncol(row_name) + 1) : (ncol(row_name) + ncol(col_name))] <- col_name
        table_m_index[temp,1] <- k1 
        #table_m_index[temp,2] <- 1:nrow(col_name)
        table_m[temp, ncol(row_name) + 1, n_id] <- round(means[k1,n_id], digits = 4)
        table_m[temp, ncol(row_name) + 2, n_id] <- pmin(round(quantile_025[k1,n_id], digits = 4), round(quantile_975[k1,n_id], digits = 4))
        table_m[temp, ncol(row_name) + 3, n_id] <- pmax(round(quantile_025[k1,n_id], digits = 4), round(quantile_975[k1,n_id], digits = 4))
      }
    }
    # reorder table_m column names, sort values column by column, keep the order of table_m_index
    table_mediator <- list(table_m = table_m, index_name = row_name, index = table_m_index, samples = est_samples)
    return(table_mediator )
    
  }