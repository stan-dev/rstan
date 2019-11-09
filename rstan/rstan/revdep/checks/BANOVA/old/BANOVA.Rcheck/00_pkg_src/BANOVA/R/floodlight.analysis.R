###
# This version, the floodlight analysis hasn't been tested for two cases: (level1_num && level1_fac) (level2_num && level1_fac)
###
# find the effects of the factor only (not including the numeric variable), using cal.mediation.effects
# find the effects of the interaction between the numeric variable and the factor, using cal.effects (specific to the numeric var), and possible modifications of effect.matrix.mediator 
#
floodlight.analysis <- 
function (sol,
          numeric_name, 
          factor_name, 
          samples_l2_param, X, Z, data = NULL, dataX = NULL, dataZ = NULL, flood_values = list()){
  
  numeric_name <- trimws(numeric_name)
  factor_name <- trimws(factor_name)
  X_names = colnames(X) # names after effect coding
  X_varNames = attr(X, "varNames") # names before coding
  Z_names = colnames(Z) # names after effect coding
  Z_varNames = attr(Z, 'varNames') # names before coding
  X_assign = attr(X, "assign")
  X_classes = attr(X, "dataClasses")
  Z_assign = attr(Z, "assign")
  Z_classes = attr(Z, "dataClasses")
  # convert samples_l2_param to est_matrix (row: X_names, column: Z_names)
  n_sample <- nrow(samples_l2_param)
  num_l1 <- length(X_assign)
  num_l2 <- length(Z_assign)
  est_matrix <- array(0 , dim = c(num_l1, num_l2, n_sample), dimnames = list(trimws(X_names), trimws(Z_names), NULL))
  for (i in 1:num_l1){
    for (j in 1:n_sample)
      est_matrix[i,,j] <- samples_l2_param[j,((i-1)*num_l2+1):((i-1)*num_l2+num_l2)]
  }
  
  # find the levels (level 1 or 2) for the numeric and categorical variables
  level1_num = numeric_name %in% X_varNames
  level2_num = numeric_name %in% Z_varNames
  level1_fac = factor_name %in% X_varNames
  level2_fac = factor_name %in% Z_varNames
  if (level1_num && level2_num) stop("The numerical variable can't appear in both levels." )
  if (level1_fac && level2_fac) stop("The categorical variable can't appear in both levels." )
  if (!level1_num && !level2_num) stop("The numerical variable isn't included in the model." )
  if (!level1_fac && !level2_fac) stop("The categorical variable isn't included in the model." )
  
  ## find the range of numeric values
  if (!is.null(data)){
    num_value <- data[,which(colnames(data) == numeric_name)]
  }else{
    if(level1_num){
      index_num_X <- which(colnames(dataX[[1]]) == numeric_name)
      num_value <- c()
      for (i in 1:length(dataX)){
        num_value <- c(num_value, dataX[[i]][, index_num_X])
      }
    }else{
      index_num_Z <- which(colnames(dataZ) == numeric_name)
      num_value <- dataZ[, index_num_Z]
    }
  }
  num_range <- range(num_value)
  
  if (level1_num && level1_fac){
    if (!((X_classes[numeric_name] == "numeric" || X_classes[numeric_name] == "integer") && X_classes[factor_name] == "factor"))
      stop('Variable types are not correct.')
    # generate var names from the factor variable (e.g. name1, name2, ..)
    ## look for the name (index) in the X_varNames
    ## look for the index in X_assign
    ## find the corresponding names in X_names
    X_names_factor <- X_names[X_assign == (which(X_varNames == factor_name) - 1)] # since assign starts from 0
    
    # generate var names from the interaction between the numeric varibale and the factor (e.g. num_name : factor1, num_name : factor2,..)
    ## use the attr(,"interactions") to find the corresponding (match the name of the two variables) index  
    X_interactions <- attr(X, "interactions_num")
    X_interactions_index <- attr(X, "interactions_numeric_index")
    inter_found_flag = F
    factor_name_map_inter_name <- list()
    if (length(X_interactions) > 0){
      for (i in 1:length(X_interactions)){
        tmp_names <- attr(X_interactions[[i]], 'names')
        if (length(tmp_names) == 2 && (numeric_name %in% tmp_names) && (factor_name %in% tmp_names)){
          inter_index <- X_interactions_index[i]
          assign_index <- which(attr(X, "assign") == inter_index)
          inter_names <- X_names[assign_index]
          split_inter_names <- strsplit(inter_names, ':')
          for (j in 1:length(split_inter_names)){
            split_inter_names[[j]] <- trimws(split_inter_names[[j]])
            # create a map, factor name : interaction name
            for (name_fac_inter in split_inter_names[[j]]){
              if (name_fac_inter %in% X_names_factor)
                factor_name_map_inter_name[[name_fac_inter]] <- inter_names[j]
            }
          }
          inter_found_flag <- T
        }
      }
    }
    if (!inter_found_flag){
      stop('No specified interaction found in the model!')
    }
    flood_eff <- cal.flood.effects(sol, est_matrix, n_sample, factor_name, numeric_name, flood_values = flood_values)

    # for each level of the factor (names found before)
    ## find these factors which is corresponding to the column of est_matrix
    ## for the first row (the intercept) of est_matrix, calculate the floodlight 
    floodlight_samples <- data.frame(array(0, dim = c(n_sample, length(X_names_factor)), dimnames = list(NULL, X_names_factor)))
    for (name_fac in X_names_factor)
      for (n_s in 1:n_sample){
        floodlight_samples[[name_fac]][n_s] = - est_matrix[name_fac, 1, n_s]/est_matrix[factor_name_map_inter_name[[name_fac]], 1, n_s]
      }
    res <- array(0, dim = c(length(X_names_factor), 3), dimnames = list(X_names_factor, c('mean', '2.5%', '97.5%')))
    for (name_fac in X_names_factor){
      res[name_fac, 1] <- mean(floodlight_samples[[name_fac]])
      res[name_fac, 2:3] <- quantile(floodlight_samples[[name_fac]], c(0.025, 0.975))
    }
    rownames(res) <-  paste(numeric_name, rownames(res), sep = ":")
    
  }
  
  if (level2_num && level1_fac){
    if (!((Z_classes[numeric_name] == "numeric" || Z_classes[numeric_name] == "integer") && X_classes[factor_name] == "factor"))
      stop('Variable types are not correct.')
    # generate var names from the factor variable (e.g. name1, name2, ..)
    ## look for the name (index) in X_varNames, Z_varNames
    ## look for the index in X_assign, Z_assign
    ## find the corresponding names in X_names, Z_names
    Z_names_num <- Z_names[Z_assign == (which(Z_varNames == numeric_name) - 1)] # since assign starts from 0
    X_names_factor <- X_names[X_assign == (which(X_varNames == factor_name) - 1)] # since assign starts from 0
    
    # for each level of the factor (names found before)
    ## find these factors which is corresponding to the column of est_matrix
    ## for the first row (the intercept) of est_matrix, calculate the floodlight 
    floodlight_samples <- data.frame(array(0, dim = c(n_sample, length(X_names_factor)), dimnames = list(NULL, X_names_factor)))
    #factor_eff_samples <- data.frame(array(0, dim = c(n_sample, length(X_names_factor)), dimnames = list(NULL, X_names_factor)))
    #numeric_eff_sample <- data.frame(array(0, dim = c(n_sample, length(Z_names_num)), dimnames = list(NULL, Z_names_num)))
    for (name_fac in X_names_factor){
      for (n_s in 1:n_sample){
        floodlight_samples[[name_fac]][n_s] = - est_matrix[name_fac, 1, n_s]/est_matrix[name_fac, Z_names_num, n_s]
      }
    }
    flood_eff <- cal.flood.effects(sol, est_matrix, n_sample, factor_name, numeric_name, flood_values = flood_values)
    res <- array(0, dim = c(length(X_names_factor), 3), dimnames = list(X_names_factor, c('mean', '2.5%', '97.5%')))
    for (name_fac in X_names_factor){
      res[name_fac, 1] <- mean(floodlight_samples[[name_fac]])
      res[name_fac, 2:3] <- quantile(floodlight_samples[[name_fac]], c(0.025, 0.975))
    }
    rownames(res) <-  paste(numeric_name, rownames(res), sep = ":")
  }
  
  if (level1_num && level2_fac){
    warning('This type of analysis is under development for a more general case. It might not be stable now.', call. = FALSE)
    if (!((X_classes[numeric_name] == "numeric" || X_classes[numeric_name] == "integer") && Z_classes[factor_name] == "factor"))
      stop('Variable types are not correct.')
    # generate var names from the factor variable (e.g. name1, name2, ..)
    ## look for the name (index) in X_varNames, Z_varNames
    ## look for the index in X_assign, Z_assign
    ## find the corresponding names in X_names, Z_names
    X_names_num <- X_names[X_assign == (which(X_varNames == numeric_name) - 1)] # since assign starts from 0
    Z_names_factor <- Z_names[Z_assign == (which(Z_varNames == factor_name) - 1)] # since assign starts from 0
    
    # for each level of the factor (names found before)
    ## find these factors which is corresponding to the column of est_matrix
    ## for the first row (the intercept) of est_matrix, calculate the floodlight 
    floodlight_samples <- data.frame(array(0, dim = c(n_sample, length(Z_names_factor)), dimnames = list(NULL, Z_names_factor)))
    for (name_fac in Z_names_factor)
      for (n_s in 1:n_sample){
        floodlight_samples[[name_fac]][n_s] = - est_matrix[1, name_fac, n_s]/est_matrix[X_names_num, name_fac, n_s]
      }
    
    flood_eff <- cal.flood.effects(sol, est_matrix, n_sample, factor_name, numeric_name, flood_values = flood_values)
    
    res <- array(0, dim = c(length(Z_names_factor), 3), dimnames = list(Z_names_factor, c('mean', '2.5%', '97.5%')))
    for (name_fac in Z_names_factor){
      res[name_fac, 1] <- mean(floodlight_samples[[name_fac]])
      res[name_fac, 2:3] <- quantile(floodlight_samples[[name_fac]], c(0.025, 0.975))
    }
    rownames(res) <-  paste(numeric_name, rownames(res), sep = ":")
  }
  
  if (level2_num && level2_fac){
    warning('This type of analysis is under development for a more general case. It might not be stable now.', call. = FALSE)
    if (!((Z_classes[numeric_name] == "numeric" || Z_classes[numeric_name] == "integer") && Z_classes[factor_name] == "factor"))
      stop('Variable types are not correct.')
    # generate var names from the factor variable (e.g. name1, name2, ..)
    ## look for the name (index) in the Z_varNames
    ## look for the index in Z_assign
    ## find the corresponding names in Z_names
    Z_names_factor <- Z_names[Z_assign == (which(Z_varNames == factor_name) - 1)] # since assign starts from 0
    
    # generate var names from the interaction between the numeric varibale and the factor (e.g. num_name : factor1, num_name : factor2,..)
    ## use the attr(,"interactions") to find the corresponding (match the name of the two variables) index  
    Z_interactions <- attr(Z, "interactions_num")
    Z_interactions_index <- attr(Z, "interactions_numeric_index")
    inter_found_flag = F
    factor_name_map_inter_name <- list()
    if (length(Z_interactions) > 0){
      for (i in 1:length(Z_interactions)){
        tmp_names <- attr(Z_interactions[[i]], 'names')
        if (length(tmp_names) == 2 && (numeric_name %in% tmp_names) && (factor_name %in% tmp_names)){
          inter_index <- Z_interactions_index[i]
          assign_index <- which(attr(Z, "assign") == inter_index)
          inter_names <- Z_names[assign_index]
          split_inter_names <- strsplit(inter_names, ':')
          for (j in 1:length(split_inter_names)){
            split_inter_names[[j]] <- trimws(split_inter_names[[j]])
            # create a map, factor name : interaction name
            for (name_fac_inter in split_inter_names[[j]]){
              if (name_fac_inter %in% Z_names_factor)
                factor_name_map_inter_name[[name_fac_inter]] <- inter_names[j]
            }
          }
          inter_found_flag <- T
        }
      }
    }
    if (!inter_found_flag){
      stop('No specified interaction found in the model!')
    }
    
    flood_eff <- cal.flood.effects(sol, est_matrix, n_sample, factor_name, numeric_name, flood_values = flood_values)
    
    # for each level of the factor (names found before)
    ## find these factors which is corresponding to the column of est_matrix
    ## for the first row (the intercept) of est_matrix, calculate the floodlight 
    floodlight_samples <- data.frame(array(0, dim = c(n_sample, length(Z_names_factor)), dimnames = list(NULL, Z_names_factor)))
    for (name_fac in Z_names_factor)
      for (n_s in 1:n_sample){
        floodlight_samples[[name_fac]][n_s] = - est_matrix[1, name_fac, n_s]/est_matrix[1, factor_name_map_inter_name[[name_fac]], n_s]
      }
    res <- array(0, dim = c(length(Z_names_factor), 3), dimnames = list(Z_names_factor, c('mean', '2.5%', '97.5%')))
    for (name_fac in Z_names_factor){
      res[name_fac, 1] <- mean(floodlight_samples[[name_fac]])
      res[name_fac, 2:3] <- quantile(floodlight_samples[[name_fac]], c(0.025, 0.975))
    }
    rownames(res) <-  paste(numeric_name, rownames(res), sep = ":")
  }
  return(list(sol = flood_eff, num_range = num_range))
}