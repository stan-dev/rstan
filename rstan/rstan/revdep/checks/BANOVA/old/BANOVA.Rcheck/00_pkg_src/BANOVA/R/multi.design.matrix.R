multi.design.matrix <-
function(l1_formula = 'NA', l2_formula = 'NA', dataX, dataZ = NULL, id, contrast = NULL){
  tmp_contrasts <- getOption("contrasts")
  #check formats
  if (l1_formula == 'NA') 
    stop("Formula in level 1 is missing! In the case of no level 1 factor, please use 'y~1'")
  
  if (l2_formula == 'NA'){
    if (class(dataX) != 'list') stop('Xsamples must be a list of features data.')
    temp1 <- nrow(dataX[[1]])
    temp2 <- ncol(dataX[[1]])
    for (i in 1:length(dataX))
      if (nrow(dataX[[i]]) != temp1 || ncol(dataX[[i]]) != temp2)
        stop('dataX dimensions mismatch! Each element of the list should have the same dimension.')
    
    n_ob <- length(dataX)
    first_choice_set <- dataX[[1]] # used to make all attributes
    # check if intercepts are included
    mf1_temp <- model.frame(formula = l1_formula, data = first_choice_set)
    intercept_bool <- attr(attr(mf1_temp,'terms'),'intercept')
    if (!intercept_bool) stop('Intercept in level 1 must be included!')
    n_choice <- nrow(first_choice_set)
    if (n_choice < 2) stop('Data error, number of choices less than 2. Please change models.')
    n_features <- ncol(first_choice_set)
    col_names <- colnames(first_choice_set)
    
    intercept_names <- ''
    for (i in 2:n_choice)
      intercept_names[i-1] <- paste('choice',i,'_Intercept', sep = "")
    #col_names <- c(intercept_names,col_names)
    
    dataClasses <- ''
    for (i in 1:ncol(first_choice_set)){
      dataClasses[i] <- class(first_choice_set[,i])
    }
    
    X_full <- list()
    X_original_choice <- list()
    # build a big matrix for each choice corresponds to all observations
    for (i in 1:n_choice){
      matrix_temp <- matrix(0, nrow = n_ob, ncol = n_features)
      for (j in 1:n_ob){
        temp <- as.matrix(dataX[[j]])
        matrix_temp[j,] <- temp[i,]
      }
      colnames(matrix_temp) <- col_names
      matrixtemp <- data.frame(matrix_temp)
      for (j in 1:n_features){
        if (dataClasses[j] == 'numeric'){
          matrixtemp[,j] <- as.numeric(as.character(matrixtemp[,j])) # R somehow automatically convert everything to factors
        }else if(dataClasses[j] == 'integer'){
          matrixtemp[,j] <- as.integer(as.character(matrixtemp[,j]))
        }else if(dataClasses[j] == 'factor'){
          matrixtemp[,j] <- as.factor(as.character(matrixtemp[,j]))
        }else 
          stop('Bad class object in dataX')
      }
      
      X_original_choice[[i]] <- matrixtemp
      
      intercept_matrix_values <- matrix(0, nrow = nrow(matrixtemp), ncol = n_choice - 1)
      if (i > 1)
        intercept_matrix_values[, i-1] <- 1
      colnames(intercept_matrix_values) <- intercept_names
      X_new_values <- cbind(intercept_matrix_values, matrixtemp) # for extracting values of variables
      
      # build full design matrix of X
      options(contrasts = rep("contr.sum",2)) # use effect coding for both unordered and ordered factors
      options(na.action='na.fail') # An error will occur if features data contains NA's.
      #### 1.1.2
      matrixtemp <- assign_contrast(matrixtemp, contrast)
      ####
      mf1 <- model.frame(formula = l1_formula, data = matrixtemp)
      #y <- model.response(mf1)
      X <- model.matrix(attr(mf1,'terms'), data = mf1)
      if (dim(X)[2] == 0) stop('No variables in level 1 model! Please add an intercept at least.')
      
      original_assign <- attr(X,'assign')
      if (attr(attr(mf1,'terms'),'intercept') != 0){
        if (ncol(X) == 2){
          tempname <- colnames(X)[2]
          X <- matrix(X[, -1],ncol = 1) # exclude the intercept because of the identification problem, this removes the attribute assign, if 2 columns before, then it becomes 1 row!
          colnames(X) <- tempname        
        }else
          X <- X[, -1]
        attr(X,'assign') <- original_assign[-1]
      }
      
      attr(X,'dataClasses') <-  attr(attr(mf1,'terms'),'dataClasses')
      numeric_index <- which(attr(X,'dataClasses') == 'numeric' | attr(X,'dataClasses') == 'integer')
      tmp <- get.interactions(mf1)
      numeric_index <- c(numeric_index, tmp$numeric_index) # add interactions including at least one numeric variable
      numeric_index_in_X <- array(NA, dim = c(1,length(numeric_index)))
      if (length(numeric_index) > 0){
        for (j in 1:length(numeric_index)){
          numeric_index_in_X[j] <- which(attr(X,'assign') == numeric_index[j])
        }
        X[,numeric_index_in_X] <- mean.center(X[,numeric_index_in_X]) # mean center numeric variables (covariates) in X
        warning("level 1 numeric variables have been mean centered.\n", call. = F, immediate. = T)
      }
      
      # then col bind the intercept matrix
      intercept_matrix <- matrix(0, nrow = n_ob, ncol = n_choice - 1)
      if (i > 1)
        intercept_matrix[, i-1] <- 1
      colnames(intercept_matrix) <- intercept_names
      X_new <- cbind(intercept_matrix, X)
      attr(X_new,'numeric_index') <- c(1:(n_choice - 1), numeric_index_in_X + n_choice - 1)
      attr(X_new,'assign') <- c(1:(n_choice - 1), attr(X,'assign') + n_choice - 1)
      attr(X_new,'varNames') <- c(intercept_names, attr(attr(mf1,'terms'),'term.labels'))
      #else
      #  stop('Intercept in level 1 must be included!')
      temp <- matrix(c(rep('numeric', n_choice - 1), attr(attr(mf1,'terms'),'dataClasses')), nrow = 1)
      colnames(temp) <- c(intercept_names, attr(attr(attr(mf1,'terms'),'dataClasses'),'names'))
      attr(X_new,'dataClasses') <- matrix(c(rep('numeric', n_choice - 1), attr(attr(mf1,'terms'),'dataClasses')), nrow = 1)
      attr(attr(X_new,'dataClasses'), 'names') <- c(intercept_names,attr(attr(attr(mf1,'terms'),'dataClasses'),'names'))
      attr(X_new,'varValues') <- get.values(data.frame(X_new_values), mf1, intercept_names = intercept_names)
      temp <- get.interactions(mf1, n_choice - 1)
      attr(X_new,'interactions') <- temp$results
      attr(X_new,'interactions_index') <- temp$index
      attr(X_new,'interactions_num') <- temp$results_num
      attr(X_new,'interactions_numeric_index') <- temp$numeric_index
      X_full[[i]] <- X_new
    }
    # check dimensions of each matrix for each choice
    n_full_feature <- ncol(X_full[[1]])
    for (i in 2:n_choice)
      if(ncol(X_full[[i]]) != n_full_feature) stop('Data is not balanced! for each choice alternative, the categorical features should have the same number of levels.')
    Z = NULL
    Z_full = NULL
  }else{
    if (class(dataX) != 'list') stop('Xsamples must be a list of features data.')
    if (length(dataX) != nrow(dataZ)) stop('X, Z samples dimension mismatch.')
    temp1 <- nrow(dataX[[1]])
    temp2 <- ncol(dataX[[1]])
    for (i in 1:length(dataX))
      if (nrow(dataX[[i]]) != temp1 || ncol(dataX[[i]]) != temp2)
        stop('dataX dimensions mismatch! Each element of the list should have the same dimension.')
      
    n_ob <- length(dataX)
    first_choice_set <- dataX[[1]] # used to make all attributes
    # check if intercepts are included
    mf1_temp <- model.frame(formula = l1_formula, data = first_choice_set)
    intercept_bool <- attr(attr(mf1_temp,'terms'),'intercept')
    if (!intercept_bool) stop('Intercept in level 1 must be included!')
    n_choice <- nrow(first_choice_set)
    if (n_choice < 2) stop('Data error, number of choices less than 2. Please change models.')
    n_features <- ncol(first_choice_set)
    col_names <- colnames(first_choice_set)
    
    intercept_names <- ''
    for (i in 2:n_choice)
      intercept_names[i-1] <- paste('choice',i,'_Intercept', sep = "")
    #col_names <- c(intercept_names,col_names)
    
    dataClasses <- ''
    for (i in 1:ncol(first_choice_set)){
      dataClasses[i] <- class(first_choice_set[,i])
    }
    
    X_full <- list()
    X_original_choice <- list()
    # build a big matrix for each choice corresponds to all observations
    for (i in 1:n_choice){
      matrix_temp <- matrix(0, nrow = n_ob, ncol = n_features)
      for (j in 1:n_ob){
        temp <- as.matrix(dataX[[j]])
        matrix_temp[j,] <- temp[i,]
      }
      colnames(matrix_temp) <- col_names
      matrixtemp <- data.frame(matrix_temp)
      for (j in 1:n_features){
        if (dataClasses[j] == 'numeric'){
          matrixtemp[,j] <- as.numeric(as.character(matrixtemp[,j])) # R somehow automatically convert everything to factors
        }else if(dataClasses[j] == 'integer'){
          matrixtemp[,j] <- as.integer(as.character(matrixtemp[,j]))
        }else if(dataClasses[j] == 'factor'){
          matrixtemp[,j] <- as.factor(as.character(matrixtemp[,j]))
        }else 
          stop('Bad class object in dataX')
      }
      
      X_original_choice[[i]] <- matrixtemp
      
      intercept_matrix_values <- matrix(0, nrow = nrow(matrixtemp), ncol = n_choice - 1)
      if (i > 1)
        intercept_matrix_values[, i-1] <- 1
      colnames(intercept_matrix_values) <- intercept_names
      X_new_values <- cbind(intercept_matrix_values, matrixtemp) # for extracting values of variables
      
      # build full design matrix of X and Z
      options(contrasts = rep("contr.sum",2)) # use effect coding for both unordered and ordered factors
      options(na.action='na.fail') # An error will occur if features data contains NA's.
      #### 1.1.2
      matrixtemp <- assign_contrast(matrixtemp, contrast)
      ####
      mf1 <- model.frame(formula = l1_formula, data = matrixtemp)
      #y <- model.response(mf1)
      X <- model.matrix(attr(mf1,'terms'), data = mf1)
      if (dim(X)[2] == 0) stop('No variables in level 1 model! Please add an intercept at least.')
    
      original_assign <- attr(X,'assign')
      if (attr(attr(mf1,'terms'),'intercept') != 0){
        if (ncol(X) == 2){
          tempname <- colnames(X)[2]
          X <- matrix(X[, -1],ncol = 1) # exclude the intercept because of the identification problem, this removes the attribute assign, if 2 columns before, then it becomes 1 row!
          colnames(X) <- tempname        
        }else
          X <- X[, -1]
        attr(X,'assign') <- original_assign[-1]
      }
    
      attr(X,'dataClasses') <-  attr(attr(mf1,'terms'),'dataClasses')
      numeric_index <- which(attr(X,'dataClasses') == 'numeric' | attr(X,'dataClasses') == 'integer')
      tmp <- get.interactions(mf1)
      numeric_index <- c(numeric_index, tmp$numeric_index) # add interactions including at least one numeric variable
      numeric_index_in_X <- array(NA, dim = c(1,length(numeric_index)))
      if (length(numeric_index) > 0){
        for (j in 1:length(numeric_index)){
          numeric_index_in_X[j] <- which(attr(X,'assign') == numeric_index[j])
        }
        X[,numeric_index_in_X] <- mean.center(X[,numeric_index_in_X]) # mean center numeric variables (covariates) in X
        warning("level 1 numeric variables have been mean centered.\n", call. = F, immediate. = T)
      }
      
      # then col bind the intercept matrix
      intercept_matrix <- matrix(0, nrow = n_ob, ncol = n_choice - 1)
      if (i > 1)
        intercept_matrix[, i-1] <- 1
      colnames(intercept_matrix) <- intercept_names
      X_new <- cbind(intercept_matrix, X)
      attr(X_new,'numeric_index') <- c(1:(n_choice - 1), numeric_index_in_X + n_choice - 1)
      attr(X_new,'assign') <- c(1:(n_choice - 1), attr(X,'assign') + n_choice - 1)
      attr(X_new,'varNames') <- c(intercept_names, attr(attr(mf1,'terms'),'term.labels'))
      #else
      #  stop('Intercept in level 1 must be included!')
      temp <- matrix(c(rep('numeric', n_choice - 1), attr(attr(mf1,'terms'),'dataClasses')), nrow = 1)
      colnames(temp) <- c(intercept_names, attr(attr(attr(mf1,'terms'),'dataClasses'),'names'))
      attr(X_new,'dataClasses') <- matrix(c(rep('numeric', n_choice - 1), attr(attr(mf1,'terms'),'dataClasses')), nrow = 1)
      attr(attr(X_new,'dataClasses'), 'names') <- c(intercept_names,attr(attr(attr(mf1,'terms'),'dataClasses'),'names'))
      attr(X_new,'varValues') <- get.values(data.frame(X_new_values), mf1, intercept_names = intercept_names)
      temp <- get.interactions(mf1, n_choice - 1)
      attr(X_new,'interactions') <- temp$results
      attr(X_new,'interactions_index') <- temp$index
      attr(X_new,'interactions_num') <- temp$results_num
      attr(X_new,'interactions_numeric_index') <- temp$numeric_index
      X_full[[i]] <- X_new
    }
    # check dimensions of each matrix for each choice
    n_full_feature <- ncol(X_full[[1]])
    for (i in 2:n_choice)
      if(ncol(X_full[[i]]) != n_full_feature) stop('Data is not balanced! for each choice alternative, the categorical features should have the same number of levels.')
    
    #### 1.1.2
    dataZ <- assign_contrast(dataZ, contrast)
    ####
    mf2 <- model.frame(formula = l2_formula, data = dataZ)
    Z_full <- model.matrix(attr(mf2,'terms'), data = mf2)
    if (dim(Z_full)[2] == 0) stop('No variables in level 2 model! Please add an intercept at least.')
    if (attr(attr(mf2,'terms'),'response') == 1) stop("level 2 model shouldn't include the response! Start with the '~'.")
    #convert Z from long format to short format, using the first row of each id
    num_id <- length(unique(id))
    id_index <- array(NA,dim = c(num_id,1))
    for (i in 1:num_id){
      index_row <- which(id == i)
      # check if it is a good between-subject factor
      if (max(index_row) > nrow(Z_full)) stop('ID indicator error!')
      for (j in 1:ncol(Z_full))
        if (max(Z_full[index_row, j]) != min(Z_full[index_row, j])) stop('Bad between-subject factor!')
      id_index[i] <- index_row[1]
    }
    Z <- as.matrix(Z_full[id_index,])
    colnames(Z) <- colnames(Z_full)
    if (attr(attr(mf2,'terms'),'intercept') != 0) # for tables in outputs
      attr(Z,'varNames') <- c('(Intercept)',attr(attr(mf2,'terms'),'term.labels'))
    else
      stop('Intercept in level 2 must be included!')
    attr(Z,'dataClasses') <- attr(attr(mf2,'terms'),'dataClasses')
    attr(Z,'varValues') <- get.values(dataZ, mf2, id_index)
    temp <- get.interactions(mf2)
    attr(Z,'interactions') <- temp$results
    attr(Z,'interactions_index') <- temp$index
    attr(Z,'interactions_num') <- temp$results_num
    attr(Z,'interactions_numeric_index') <- temp$numeric_index
    attr(Z,"assign") <- attr(Z_full,"assign")
    attr(Z,"contrasts") <- attr(Z_full,"contrasts")
    # find index of numeric variables (covariates) in Z
    numeric_index <- which(attr(Z,'dataClasses') == 'numeric' | attr(Z,'dataClasses') == 'integer')
    numeric_index <- c(numeric_index, temp$numeric_index) # add interactions including at least one numeric variable
    numeric_index_in_Z <- array(NA, dim = c(1,length(numeric_index)))
    if (length(numeric_index) > 0){
      numeric_index_in_Z <- array(NA, dim = c(1,length(numeric_index)))
      for (i in 1:length(numeric_index)){
        numeric_index_in_Z[i] <- which(attr(Z,'assign') == numeric_index[i])
      }
      Z[,numeric_index_in_Z] <- mean.center(Z[,numeric_index_in_Z])
      warning("level 2 numeric variables have been mean centered.\n", call. = F, immediate. = T)
    }
    attr(Z,'numeric_index') <- numeric_index_in_Z
    attr(Z_full,'numeric_index') <- numeric_index_in_Z
  }
  options(contrasts = tmp_contrasts)
  return(list(X_full = X_full, X_original_choice = X_original_choice, Z = Z, Z_full = Z_full))
}
