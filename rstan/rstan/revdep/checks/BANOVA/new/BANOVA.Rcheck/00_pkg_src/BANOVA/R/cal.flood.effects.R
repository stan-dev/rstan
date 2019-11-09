# this function outputs effects of the numeric variable (sum of all coefficients related to 
# the interaction of the numeric variable and the factor for each combination of moderators)in the model sol_1
# calculate the interaction effects between the numeric variable and the factor at different levels of moderators 

cal.flood.effects <-
  function (sol_1, est_matrix, n_sample, factor_name, numeric_name, flood_values = list()){
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
    factor_in_l1 <- factor_name %in% rownames(model1_level1_var_matrix)
    factor_in_l2 <- factor_name %in% rownames(model1_level2_var_matrix)
    numeric_in_l1 <- numeric_name %in% rownames(model1_level1_var_matrix)
    numeric_in_l2 <- numeric_name %in% rownames(model1_level2_var_matrix)
    if (!factor_in_l1 & !factor_in_l2) stop(factor_name," is not included in the model!")
    if (!numeric_in_l1 & !numeric_in_l2) stop(numeric_name," is not included in the model!")
    # mediator is at level 1
    if (factor_in_l1 & numeric_in_l2 ){
      factor_assign <- which (rownames(model1_level1_var_matrix) == factor_name)
      #numeric_assign <- which (rownames(model1_level2_var_matrix) == numeric_name)
      interaction_list <- attr(sol_1$dMatrice$X, "interactions_num") # if the mediator is numeric 
      interaction_list_index <- attr(sol_1$dMatrice$X, "interactions_numeric_index")
      factor_interaction_list <- list()
      factor_interaction_list_index <- array()
      j = 1
      if (length(interaction_list) > 0){
        for (i in 1:length(interaction_list)){
          if (factor_assign %in% interaction_list[[i]]){
            factor_interaction_list[[j]] = interaction_list[[i]]
            factor_interaction_list_index[j] = interaction_list_index[i]
            j = j + 1
          }
        }
      }
      interaction_list <- attr(sol_1$dMatrice$X, "interactions") # if the mediator is not numeric 
      interaction_list_index <- attr(sol_1$dMatrice$X, "interactions_index")
      if (length(interaction_list) > 0){
        for (i in 1:length(interaction_list)){
          if (factor_assign %in% interaction_list[[i]]){
            factor_interaction_list[[j]] = interaction_list[[i]]
            factor_interaction_list_index[j] = interaction_list_index[i]
            j = j + 1
          }
        }
      }
      ################################################
      # find the effects of the factor and the numeric
      ################################################
      # for each interaction in factor_interaction_list create a moderated table (exclude main effects)
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
        num_l2_matrix <- effect.matrix.mediator(interaction_factors = l2_values,
                                                matrix_formula=formula(attr(sol_1$mf2, 'terms')),
                                                mediator=numeric_name, 
                                                flood_values = flood_values, contrast = sol_1$contrast)
        l2_matrix <- effect.matrix.mediator(interaction_factors = l2_values, 
                                            matrix_formula=formula(attr(sol_1$mf2, 'terms')), 
                                            xvar=numeric_name, intercept_include = TRUE, 
                                            flood_values = flood_values, contrast = sol_1$contrast)
      }
      if (length(factor_interaction_list) > 0){
        l1_values <- attr(sol_1$dMatrice$X, 'varValues')
        factor_interaction_effect_matrix <- list()
        for (i in 1:length(factor_interaction_list)){
          # l1 matrix
          # TO check:
          # y var is also included in l1_values, interaction_list has considered this
          factor_interaction_effect_matrix[[i]] <- effect.matrix.mediator(interaction_factors = l1_values[factor_interaction_list[[i]]], 
                                                                          mediator=factor_name, xvar=numeric_name, 
                                                                          flood_values = flood_values, contrast = sol_1$contrast)
          est_samples <- array(0, dim = c(nrow(factor_interaction_effect_matrix[[i]]), nrow(l2_matrix), n_sample))
          num_est_samples <- array(0, dim = c(nrow(factor_interaction_effect_matrix[[i]]), nrow(num_l2_matrix), n_sample))
          for (n_s in 1:n_sample){
            est_samples[,,n_s] <- factor_interaction_effect_matrix[[i]] %*% est_matrix[colnames(factor_interaction_effect_matrix[[i]]), colnames(l2_matrix), n_s] %*% t(l2_matrix)
            num_est_samples[,,n_s] <- factor_interaction_effect_matrix[[i]] %*% est_matrix[colnames(factor_interaction_effect_matrix[[i]]), colnames(num_l2_matrix), n_s] %*% t(num_l2_matrix)
          }
          #table_m <- construct.table(est_samples, attr(factor_interaction_effect_matrix[[i]], 'levels'), attr(l2_matrix, 'levels'))
          table_m <- construct.effect.table(est_samples, attr(factor_interaction_effect_matrix[[i]], 'levels'), attr(l2_matrix, 'levels'), 
                                            num_est_samples, attr(factor_interaction_effect_matrix[[i]], 'levels'), attr(num_l2_matrix, 'levels'), l1_values[factor_assign], numeric_name)
          table_list[[table_list_index]] <- table_m
          table_list_index = table_list_index + 1
        }
        
      }else{
        l1_values <- attr(sol_1$dMatrice$X, 'varValues')
        # no interaction with the mediator in level 1, only select the mediator
        # TO check:
        # y var is also included in l1_values, interaction_list has considered this
        l1_matrix <- effect.matrix.mediator(l1_values[factor_assign], 
                                            mediator=factor_name, 
                                            flood_values = flood_values, contrast = sol_1$contrast)
        est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(l2_matrix), n_sample))
        num_est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(num_l2_matrix), n_sample))
        for (n_s in 1:n_sample){
          est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(l2_matrix), n_s] %*% t(l2_matrix)
          num_est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(num_l2_matrix), n_s] %*% t(num_l2_matrix)
        }
        #table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'))
        table_m <- construct.effect.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'), 
                                          num_est_samples, attr(l1_matrix, 'levels'), attr(num_l2_matrix, 'levels'), l1_values[factor_assign], numeric_name)
        table_list[[table_list_index]] <- table_m
        table_list_index = table_list_index + 1
      }
      
    }
    
    if (factor_in_l1 & numeric_in_l1 ){
      # find if there are moderators interacts with the factor except the numeric
      factor_assign <- which (rownames(model1_level1_var_matrix) == factor_name)
      numeric_assign <- which (rownames(model1_level1_var_matrix) == numeric_name)
      interaction_list <- attr(sol_1$dMatrice$X, "interactions_num") # if the factor is numeric 
      interaction_list_index <- attr(sol_1$dMatrice$X, "interactions_numeric_index")
      factor_interaction_list <- list()
      factor_interaction_list_index <- array()
      j = 1
      if (length(interaction_list) > 0){
        for (i in 1:length(interaction_list)){
          if (factor_assign %in% interaction_list[[i]] && !(numeric_assign %in% interaction_list[[i]])){
            factor_interaction_list[[j]] = interaction_list[[i]]
            factor_interaction_list_index[j] = interaction_list_index[i]
            j = j + 1
          }
        }
      }
      interaction_list <- attr(sol_1$dMatrice$X, "interactions") # if the factor is not numeric
      interaction_list_index <- attr(sol_1$dMatrice$X, "interactions_index")
      if (length(interaction_list) > 0){
        for (i in 1:length(interaction_list)){
          if (factor_assign %in% interaction_list[[i]]){
            factor_interaction_list[[j]] = interaction_list[[i]]
            factor_interaction_list_index[j] = interaction_list_index[i]
            j = j + 1
          }
        }
      }
      # for each interaction in factor_interaction_list create a moderated table (exclude main effects)
      # calculate l2 matrix first using all params but xvar, also using l2 formula
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
        l2_matrix <- effect.matrix.mediator(interaction_factors = l2_values, 
                                            matrix_formula=formula(attr(sol_1$mf2, 'terms')),
                                            flood_values = flood_values, contrast = sol_1$contrast)
      }
      if (length(factor_interaction_list) > 0){
        #TODO: check if numeric varible in L1
        l1_values <- attr(sol_1$dMatrice$X, 'varValues')
        factor_interaction_effect_matrix <- list()
        num_factor_interaction_effect_matrix <- list()
        # remove y
        # find the response variable
        # vars <- as.character(attr(attr(sol_1$mf1, 'terms'), 'variables'))[-1]
        # response_ind <- attr(attr(sol_1$mf1, 'terms'), 'response')
        # response_name <- vars[response_ind]
        # #remove the response value from l1_values
        # to_rm_res <- 0
        # for (i in 1:length(l1_values)){
        #   if (attr(l1_values[[i]], 'var_names') == response_name) to_rm_res <- i
        # }
        # num_l1_values <- l1_values
        # num_l1_values[[to_rm_res]] <- NULL
        for (i in 1:length(factor_interaction_list)){
          # l1 matrix should include the effects of interactions contains both the numeric and factor
          # TO check:
          # y var is also included in l1_values, interaction_list has considered this
          factor_interaction_effect_matrix[[i]] <- effect.matrix.mediator(interaction_factors = l1_values[factor_interaction_list[[i]]], 
                                                                          mediator=factor_name, 
                                                                          xvar=numeric_name, 
                                                                          flood_values = flood_values,
                                                                          contrast = sol_1$contrast)
          num_factor_interaction_effect_matrix[[i]] <- effect.matrix.mediator(interaction_factors = l1_values,
                                                                              matrix_formula=formula(attr(sol_1$mf1, 'terms')),
                                                                              mediator=factor_name, 
                                                                              xvar=numeric_name, 
                                                                              xvar_include = TRUE, 
                                                                              flood_values = flood_values,
                                                                              contrast = sol_1$contrast)
          est_samples <- array(0, dim = c(nrow(factor_interaction_effect_matrix[[i]]), nrow(l2_matrix), n_sample))
          num_est_samples <- array(0, dim = c(nrow(num_factor_interaction_effect_matrix[[i]]), nrow(l2_matrix), n_sample))
          if ("" %in% dimnames(est_matrix)[[2]]){
            for (n_s in 1:n_sample){
              est_samples[,,n_s] <- factor_interaction_effect_matrix[[i]] %*% est_matrix[colnames(factor_interaction_effect_matrix[[i]]), 1, n_s] %*% t(l2_matrix)
              num_est_samples[,,n_s] <- num_factor_interaction_effect_matrix[[i]] %*% est_matrix[colnames(num_factor_interaction_effect_matrix[[i]]), 1, n_s] %*% t(l2_matrix)
            }
          }else{
            for (n_s in 1:n_sample){
              est_samples[,,n_s] <- factor_interaction_effect_matrix[[i]] %*% est_matrix[colnames(factor_interaction_effect_matrix[[i]]), colnames(l2_matrix), n_s] %*% t(l2_matrix)
              num_est_samples[,,n_s] <- num_factor_interaction_effect_matrix[[i]] %*% est_matrix[colnames(num_factor_interaction_effect_matrix[[i]]), colnames(l2_matrix), n_s] %*% t(l2_matrix)
            }
          }
          #table_m <- construct.table(est_samples, attr(factor_interaction_effect_matrix[[i]], 'levels'), attr(l2_matrix, 'levels'))
          table_m <- construct.effect.table(est_samples, attr(factor_interaction_effect_matrix[[i]], 'levels'), attr(l2_matrix, 'levels'), 
                                            num_est_samples, attr(num_factor_interaction_effect_matrix[[i]], 'levels'), attr(l2_matrix, 'levels'), l1_values[factor_assign], numeric_name)
          table_list[[table_list_index]] <- table_m
          table_list_index = table_list_index + 1
        }
        
      }else{
        l1_values <- attr(sol_1$dMatrice$X, 'varValues')
        factor_values <- l1_values[factor_assign]
        # no interaction with the mediator in level 1, only select the mediator
        # check:
        # y var is also included in l1_values, interaction_list has considered this
        l1_matrix <- effect.matrix.mediator(l1_values[factor_assign], 
                                            mediator=factor_name, 
                                            flood_values = flood_values, contrast = sol_1$contrast)
        # find the response variable
        # vars <- as.character(attr(attr(sol_1$mf1, 'terms'), 'variables'))[-1]
        # response_ind <- attr(attr(sol_1$mf1, 'terms'), 'response')
        # response_name <- vars[response_ind]
        # #remove the response value from l1_values
        # to_rm_res <- 0
        # for (i in 1:length(l1_values)){
        #   if (attr(l1_values[[i]], 'var_names') == response_name) to_rm_res <- i
        # }
        # l1_values[[to_rm_res]] <- NULL
        num_l1_matrix <- effect.matrix.mediator(l1_values, 
                                                matrix_formula=formula(attr(sol_1$mf1, 'terms')),
                                                mediator=factor_name, 
                                                xvar = numeric_name, 
                                                xvar_include = TRUE, 
                                                flood_values = flood_values, 
                                                contrast = sol_1$contrast)
        est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(l2_matrix), n_sample))
        num_est_samples <- array(0, dim = c(nrow(num_l1_matrix), nrow(l2_matrix), n_sample))
        if ("" %in% dimnames(est_matrix)[[2]]){
          for (n_s in 1:n_sample){
            est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), 1, n_s] %*% t(l2_matrix)
            num_est_samples[,,n_s] <- num_l1_matrix %*% est_matrix[colnames(num_l1_matrix), 1, n_s] %*% t(l2_matrix)
          }
        }else{
          for (n_s in 1:n_sample){
            est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(l2_matrix), n_s] %*% t(l2_matrix)
            num_est_samples[,,n_s] <- num_l1_matrix %*% est_matrix[colnames(num_l1_matrix), colnames(l2_matrix), n_s] %*% t(l2_matrix)
          }
        }
        
        #table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'))
        table_m <- construct.effect.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'), 
                                          num_est_samples, attr(num_l1_matrix, 'levels'), attr(l2_matrix, 'levels'), factor_values, numeric_name)
        table_list[[table_list_index]] <- table_m
        table_list_index = table_list_index + 1
      }
      
    }
  
    if (factor_in_l2 & numeric_in_l1 ){
      # find if there are moderators interacts with the factor
      factor_assign <- which (rownames(model1_level2_var_matrix) == factor_name)
      #numeric_assign <- which (rownames(model1_level1_var_matrix) == numeric_name)
      interaction_list <- attr(sol_1$dMatrice$Z, "interactions_num")
      interaction_list_index <- attr(sol_1$dMatrice$Z, "interactions_numeric_index")
      factor_interaction_list <- list()
      factor_interaction_list_index <- array()
      j = 1
      if (length(interaction_list) > 0){
        for (i in 1:length(interaction_list)){
          if (factor_assign %in% interaction_list[[i]]){
            factor_interaction_list[[j]] = interaction_list[[i]]
            factor_interaction_list_index[j] = interaction_list_index[i]
            j = j + 1
          }
        }
      }
      interaction_list <- attr(sol_1$dMatrice$X, "interactions") # if the factor is not numeric
      interaction_list_index <- attr(sol_1$dMatrice$X, "interactions_index")
      if (length(interaction_list) > 0){
        for (i in 1:length(interaction_list)){
          if (factor_assign %in% interaction_list[[i]]){
            factor_interaction_list[[j]] = interaction_list[[i]]
            factor_interaction_list_index[j] = interaction_list_index[i]
            j = j + 1
          }
        }
      }
      
      l1_values <- attr(sol_1$dMatrice$X, 'varValues')
      # vars <- as.character(attr(attr(sol_1$mf1, 'terms'), 'variables'))[-1]
      # response_ind <- attr(attr(sol_1$mf1, 'terms'), 'response')
      # response_name <- vars[response_ind]
      # #remove the response value from l1_values
      # to_rm_res <- 0
      # for (i in 1:length(l1_values)){
      #   if (attr(l1_values[[i]], 'var_names') == response_name) to_rm_res <- i
      # }
      # l1_values[[to_rm_res]] <- NULL
      
      if (length(l1_values) == 0){
        # only intercept included
        l1_matrix <- model.matrix(~1)
        l1_matrix <- rbind(l1_matrix, c(1))
        attr(l1_matrix, "levels") <- l1_matrix
      }else{
        num_l1_matrix <- effect.matrix.mediator(interaction_factors = l1_values, 
                                                matrix_formula=formula(attr(sol_1$mf1, 'terms')), 
                                                mediator=numeric_name, 
                                                flood_values = flood_values, contrast = sol_1$contrast)
        l1_matrix <- effect.matrix.mediator(interaction_factors = l1_values, 
                                            matrix_formula=formula(attr(sol_1$mf1, 'terms')), 
                                            xvar=numeric_name, 
                                            intercept_include = TRUE, 
                                            flood_values = flood_values, 
                                            contrast = sol_1$contrast)
      }
      l2_values <- attr(sol_1$dMatrice$Z, 'varValues')
      if (length(l2_values) == 0){
        # only intercept included
        num_l2_matrix <- model.matrix(~1)
        num_l2_matrix <- rbind(num_l2_matrix, c(1))
        attr(num_l2_matrix, "levels") <- num_l2_matrix
        if (sol_1$single_level){
          colnames(num_l2_matrix) <- c(" ")
        }
      }else{
        num_l2_matrix <- effect.matrix.mediator(interaction_factors = l2_values, 
                                                matrix_formula=formula(attr(sol_1$mf2, 'terms')), 
                                                mediator = factor_name, 
                                                flood_values = flood_values, 
                                                contrast = sol_1$contrast)
      }
      num_est_samples <- array(0, dim = c(nrow(num_l1_matrix), nrow(num_l2_matrix), n_sample))
      for (n_s in 1:n_sample){
        num_est_samples[,,n_s] <- num_l1_matrix %*% est_matrix[colnames(num_l1_matrix), colnames(num_l2_matrix), n_s] %*% t(num_l2_matrix)
      }
      if (length(factor_interaction_list) > 0){
        #l2_values <- attr(sol_1$dMatrice$Z, 'varValues')
        factor_interaction_effect_matrix <- list()
        for (i in 1:length(factor_interaction_list)){
          # l2 matrix
          factor_interaction_effect_matrix[[i]] <- effect.matrix.mediator(interaction_factors = l2_values[factor_interaction_list[[i]]], 
                                                                          mediator=factor_name, 
                                                                          xvar=numeric_name, 
                                                                          flood_values = flood_values, 
                                                                          contrast = sol_1$contrast)
          est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(factor_interaction_effect_matrix[[i]]), n_sample))
          for (n_s in 1:n_sample){
            est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(factor_interaction_effect_matrix[[i]]), n_s] %*% t(factor_interaction_effect_matrix[[i]])
          }
          #TODO use a list to store
          #table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(mediator_interaction_effect_matrix[[i]], 'levels'))
          table_m <- construct.effect.table(est_samples, attr(l1_matrix, 'levels'), attr(factor_interaction_effect_matrix[[i]], 'levels'), 
                                            num_est_samples, attr(num_l1_matrix, 'levels'), attr(num_l2_matrix, 'levels'), l2_values[factor_assign], numeric_name)
          table_list[[table_list_index]] <- table_m
          table_list_index = table_list_index + 1
        }
      }else{
        l2_values <- attr(sol_1$dMatrice$Z, 'varValues')
        # no interaction with the mediator in level 2, only select the mediator
        l2_matrix <- effect.matrix.mediator(l2_values[factor_assign], 
                                            mediator= factor_name, 
                                            flood_values = flood_values, contrast = sol_1$contrast)
        est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(l2_matrix), n_sample))
        for (n_s in 1:n_sample){
          est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(l2_matrix), n_s] %*% t(l2_matrix)
        }
        #table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'))
        table_m <- construct.effect.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'), 
                                          num_est_samples, attr(num_l1_matrix, 'levels'), attr(num_l2_matrix, 'levels'), l2_values[factor_assign], numeric_name)
        table_list[[table_list_index]] <- table_m
        table_list_index = table_list_index + 1
      }
      
      # est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(l2_matrix), n_sample))
      # for (n_s in 1:n_sample){
      #   est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(l2_matrix), n_s] %*% t(l2_matrix)
      # }
      # table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'))
      # table_list[[table_list_index]] <- table_m
      # table_list_index = table_list_index + 1
  
    }
    
    if (factor_in_l2 & numeric_in_l2 ){
      # find if there are moderators interacts with the factor except the numeric
      # find if there are moderators interacts with the factor
      factor_assign <- which (rownames(model1_level2_var_matrix) == factor_name)
      numeric_assign <- which (rownames(model1_level2_var_matrix) == numeric_name)
      interaction_list <- attr(sol_1$dMatrice$Z, "interactions_num")
      interaction_list_index <- attr(sol_1$dMatrice$Z, "interactions_numeric_index")
      factor_interaction_list <- list()
      factor_interaction_list_index <- array()
      j = 1
      if (length(interaction_list) > 0){
        for (i in 1:length(interaction_list)){
          if (factor_assign %in% interaction_list[[i]] && !(numeric_assign %in% interaction_list[[i]])){
            factor_interaction_list[[j]] = interaction_list[[i]]
            factor_interaction_list_index[j] = interaction_list_index[i]
            j = j + 1
          }
        }
      }
      interaction_list <- attr(sol_1$dMatrice$X, "interactions") # if the factor is not numeric
      interaction_list_index <- attr(sol_1$dMatrice$X, "interactions_index")
      if (length(interaction_list) > 0){
        for (i in 1:length(interaction_list)){
          if (factor_assign %in% interaction_list[[i]]){
            factor_interaction_list[[j]] = interaction_list[[i]]
            factor_interaction_list_index[j] = interaction_list_index[i]
            j = j + 1
          }
        }
      }

      l1_values <- attr(sol_1$dMatrice$X, 'varValues')
      # exclude y var
      # find the response variable
      # vars <- as.character(attr(attr(sol_1$mf1, 'terms'), 'variables'))[-1]
      # response_ind <- attr(attr(sol_1$mf1, 'terms'), 'response')
      # response_name <- vars[response_ind]
      # #remove the response value from l1_values
      # to_rm_res <- 0
      # for (i in 1:length(l1_values)){
      #   if (attr(l1_values[[i]], 'var_names') == response_name) to_rm_res <- i
      # }
      # l1_values[[to_rm_res]] <- NULL
      
      if (length(l1_values) == 0){
        # only intercept included
        l1_matrix <- model.matrix(~1)
        l1_matrix <- rbind(l1_matrix, c(1))
        attr(l1_matrix, "levels") <- l1_matrix
      }else{
        l1_matrix <- effect.matrix.mediator(interaction_factors = l1_values, 
                                            matrix_formula=formula(attr(sol_1$mf1, 'terms')),
                                            flood_values = flood_values, contrast = sol_1$contrast)
      }
      l2_values <- attr(sol_1$dMatrice$Z, 'varValues')
      num_l2_matrix <- effect.matrix.mediator(l2_values, 
                                              matrix_formula=formula(attr(sol_1$mf2, 'terms')),
                                              mediator=factor_name, 
                                              xvar = numeric_name, 
                                              xvar_include = TRUE, 
                                              flood_values = flood_values, 
                                              contrast = sol_1$contrast)
      num_est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(num_l2_matrix), n_sample))
      for (n_s in 1:n_sample)
        num_est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(num_l2_matrix), n_s] %*% t(num_l2_matrix)
      if (length(factor_interaction_list) > 0){
        #l2_values <- attr(sol_1$dMatrice$Z, 'varValues')
        factor_interaction_effect_matrix <- list()
        for (i in 1:length(factor_interaction_list)){
          # l2 matrix
          factor_interaction_effect_matrix[[i]] <- effect.matrix.mediator(interaction_factors = l2_values[factor_interaction_list[[i]]], 
                                                                          mediator=factor_name, 
                                                                          xvar=numeric_name, 
                                                                          flood_values = flood_values, 
                                                                          contrast = sol_1$contrast)
          est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(factor_interaction_effect_matrix[[i]]), n_sample))
          for (n_s in 1:n_sample){
            est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(factor_interaction_effect_matrix[[i]]), n_s] %*% t(factor_interaction_effect_matrix[[i]])
          }
          #TODO use a list to store
          #table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(mediator_interaction_effect_matrix[[i]], 'levels'))
          table_m <- construct.effect.table(est_samples, attr(l1_matrix, 'levels'), attr(factor_interaction_effect_matrix[[i]], 'levels'), 
                                            num_est_samples, attr(l1_matrix, 'levels'), attr(num_l2_matrix, 'levels'), l2_values[factor_assign], numeric_name)
          table_list[[table_list_index]] <- table_m
          table_list_index = table_list_index + 1
        }
      }else{
        l2_values <- attr(sol_1$dMatrice$Z, 'varValues')
        # no interaction with the mediator in level 2, only select the mediator
        l2_matrix <- effect.matrix.mediator(l2_values[factor_assign], 
                                            mediator= factor_name, 
                                            flood_values = flood_values, contrast = sol_1$contrast)
        est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(l2_matrix), n_sample))
        for (n_s in 1:n_sample){
          est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(l2_matrix), n_s] %*% t(l2_matrix)
        }
        #table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'))
        table_m <- construct.effect.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'), 
                                          num_est_samples, attr(l1_matrix, 'levels'), attr(num_l2_matrix, 'levels'), l2_values[factor_assign], numeric_name)
        table_list[[table_list_index]] <- table_m
        table_list_index = table_list_index + 1
      }    
      
      # l2_values <- attr(sol_1$dMatrice$Z, 'varValues')
      # num_l2_matrix <- effect.matrix.mediator(l2_values, mediator=factor_name, xvar = numeric_name, xvar_include = TRUE)
      # num_est_samples <- array(0, dim = c(nrow(l1_matrix), nrow(num_l2_matrix), n_sample))
      # for (n_s in 1:n_sample)
      #   num_est_samples[,,n_s] <- l1_matrix %*% est_matrix[colnames(l1_matrix), colnames(num_l2_matrix), n_s] %*% t(num_l2_matrix)
      # 
      # table_m <- construct.table(est_samples, attr(l1_matrix, 'levels'), attr(l2_matrix, 'levels'))
      # table_list[[table_list_index]] <- table_m
      # table_list_index = table_list_index + 1
    }

    return(table_list)
  }

construct.effect.table <- 
  function (est_samples, row_name, col_name, num_est_samples, num_row_name, num_col_name, factor_values, numeric_name){
    n_sample <- dim(est_samples)[3]
    sample_names <- paste('sample_', 1:n_sample, sep = "")
    table_f <- generate.table.samples(est_samples, row_name, col_name)
    table_num <- generate.table.samples(num_est_samples, num_row_name, num_col_name)
    names_to_merge <- intersect(c(colnames(row_name), colnames(col_name)), c(colnames(num_row_name), colnames(num_col_name)))
    samples_tb <- merge(table_f$tb, table_num$tb, by = names_to_merge, all.x = TRUE, suffixes = c(".x",".y"))
    names_union <-  union(c(colnames(row_name), colnames(col_name)), c(colnames(num_row_name), colnames(num_col_name)))
    # calculate floodlights
    floodlight_table <- array(NA, dim = c(nrow(samples_tb), length(names_union) + 3), dimnames = list(rep("",nrow(samples_tb)), c(names_union,'mean', '2.5%', '97.5%')))
    floodlight_table[, names_union] <- as.matrix(samples_tb[, names_union])
    for (i in 1: nrow(samples_tb)){
      sample_values <- -as.numeric(as.matrix(samples_tb[i, paste(sample_names, ".x", sep = "")]))/as.numeric(as.matrix(samples_tb[i, paste(sample_names, ".y", sep = "")]))
      floodlight_table[i, 'mean'] <- round(mean(sample_values, na.rm = TRUE), digits = 4)
      floodlight_table[i, '2.5%'] <- round(quantile(sample_values, probs = 0.025, na.rm = TRUE), digits = 4)
      floodlight_table[i, '97.5%'] <- round(quantile(sample_values, probs = 0.975, na.rm = TRUE), digits = 4)
    }
    if (length(factor_values) > 1) stop("The format of the factor values error!")
    factor_name <- attr(factor_values[[1]], "var_names")
    # select only one value of the factor
    num_levels <- length(unique(factor_values[[1]]))
    if (num_levels > 2) stop("The number of levels of the factor is greater than 2. This is not supported yet.", call. = FALSE)
    if (num_levels == 1) stop("The number of levels of the factor should be equal to two.", call. = FALSE)
    floodlight_table <- subset(floodlight_table, floodlight_table[, factor_name] == as.character(factor_values[[1]][1]))
    if (numeric_name %in% colnames(floodlight_table))  floodlight_table <- floodlight_table[, -which(colnames(floodlight_table) == numeric_name), drop = FALSE]
    if (factor_name %in% colnames(floodlight_table)) floodlight_table <- floodlight_table[, -which(colnames(floodlight_table) == factor_name), drop = FALSE]
    if ("(Intercept)" %in% colnames(floodlight_table)) floodlight_table <- floodlight_table[, -which(colnames(floodlight_table) == "(Intercept)"), drop = FALSE]
    # attach a column named numeric : factor
    rownames(floodlight_table) <- rep(paste(factor_name, ":", numeric_name), nrow(floodlight_table))
    return(floodlight_table)
    
  }

generate.table.samples <- function(est_samples, row_name, col_name){
  n_sample <- dim(est_samples)[3]
  table_f <- array(NA, dim = c(nrow(row_name) * nrow(col_name), ncol(row_name) + ncol(col_name) + n_sample), 
                   dimnames = list(rep("", nrow(row_name) * nrow(col_name)), c(colnames(row_name), colnames(col_name), paste('sample_', 1:n_sample, sep = ""))))
  table_f_index <- array(NA, dim = c(nrow(row_name) * nrow(col_name), 2), 
                         dimnames = list(NULL, c('est_samples_row_index', 'est_samples_col_index')))
  for (k1 in 1:nrow(row_name)){
    temp <- ((k1-1) * nrow(col_name) + 1):((k1-1) * nrow(col_name) + nrow(col_name))
    table_f[temp, 1:ncol(row_name)] <- t(replicate(nrow(col_name), row_name[k1, ]))
    table_f[temp, (ncol(row_name) + 1) : (ncol(row_name) + ncol(col_name))] <- col_name
    table_f_index[temp,1] <- k1
    table_f_index[temp,2] <- 1:nrow(col_name)
    r_ind <- temp[1] - 1
    for (c_ind in 1:nrow(col_name)){
      for(s_ind in 1:n_sample){
        table_f[r_ind + c_ind, paste('sample_',s_ind, sep ="")] <- est_samples[k1,c_ind,s_ind]
      }
    }
  }
  return(list(tb = table_f, tb_ind <- table_f_index))
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
