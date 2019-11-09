design.matrix <-
function(l1_formula = 'NA', l2_formula = 'NA', data, id, contrast = NULL){
  tmp_contrasts <- getOption("contrasts")
  # error checking
  if (l1_formula == 'NA') 
    stop("Formula in level 1 is missing! In the case of no level 1 factor, please use 'y~1'")
  #### 1.1.2
  data <- assign_contrast(data, contrast)
  ####
  if (l2_formula == 'NA'){
    # stop("Formula in level 2 is missing! In the case of no level 1 factor, please use '~1'")
    # one level model 
    options(contrasts = rep("contr.sum",2)) # use effect coding for both unordered and ordered factors
    options(na.action='na.pass') # pass NAs, then check features but not the response(since it's Bayesian)
    mf1 <- model.frame(formula = l1_formula, data = data)
    y <- model.response(mf1)
    if (length(y) == 0) stop('The response variable is not found!')
    X <- model.matrix(attr(mf1,'terms'), data = mf1)
    if (sum(is.na(X)) > 0) stop('Missing/NA values in level 1 features!')
    if (dim(X)[2] == 0) stop('No variables in level 1 model! Please add an intercept at least.')
    if (attr(attr(mf1,'terms'),'intercept') != 0) # for tables in outputs
      attr(X,'varNames') <- c('(Intercept)',attr(attr(mf1,'terms'),'term.labels')) 
    else
      stop('Intercept in level 1 must be included!')
    attr(X,'dataClasses') <- attr(attr(mf1,'terms'),'dataClasses')[-1] # exclude y
    attr(X,'varValues') <- get.values(data, mf1)
    temp <- get.interactions(mf1)
    attr(X,'interactions') <- temp$results
    attr(X,'interactions_index') <- temp$index
    attr(X,'interactions_num') <- temp$results_num
    attr(X,'interactions_numeric_index') <- temp$numeric_index
    # find index of numeric variables (covariates) in X
    numeric_index <- which(attr(X,'dataClasses') == 'numeric' | attr(X,'dataClasses') == 'integer')
    numeric_index <- c(numeric_index, temp$numeric_index) # add interactions including at least one numeric variable
    numeric_index_in_X <- array(NA, dim = c(1,length(numeric_index)))
    if (length(numeric_index) > 0){
      for (i in 1:length(numeric_index)){
        numeric_index_in_X[i] <- which(attr(X,'assign') == numeric_index[i])
      }
      X[,numeric_index_in_X] <- mean.center(X[,numeric_index_in_X])
      warning("level 1 numeric variables have been mean centered.\n", call. = F, immediate. = T)
    }
    attr(X,'numeric_index') <- numeric_index_in_X
    Z = NULL
    Z_full = NULL
  }else{
    options(contrasts = rep("contr.sum",2)) # use effect coding for both unordered and ordered factors
    options(na.action='na.pass') # pass NAs, then check features but not the response(since it's Bayesian)
    mf1 <- model.frame(formula = l1_formula, data = data)
    y <- model.response(mf1)
    if (length(y) == 0) stop('The response variable is not found!')
    X <- model.matrix(attr(mf1,'terms'), data = mf1)
    if (sum(is.na(X)) > 0) stop('Missing/NA values in level 1 features!')
    if (dim(X)[2] == 0) stop('No variables in level 1 model! Please add an intercept at least.')
    if (attr(attr(mf1,'terms'),'intercept') != 0) # for tables in outputs
      attr(X,'varNames') <- c('(Intercept)',attr(attr(mf1,'terms'),'term.labels')) 
    else
      stop('Intercept in level 1 must be included!')
    attr(X,'dataClasses') <- attr(attr(mf1,'terms'),'dataClasses')[-1] # exclude y
    attr(X,'varValues') <- get.values(data, mf1)
    temp <- get.interactions(mf1)
    attr(X,'interactions') <- temp$results
    attr(X,'interactions_index') <- temp$index
    attr(X,'interactions_num') <- temp$results_num
    attr(X,'interactions_numeric_index') <- temp$numeric_index
    # find index of numeric variables (covariates) in X
    numeric_index <- which(attr(X,'dataClasses') == 'numeric' | attr(X,'dataClasses') == 'integer')
    numeric_index <- c(numeric_index, temp$numeric_index) # add interactions including at least one numeric variable
    numeric_index_in_X <- array(NA, dim = c(1,length(numeric_index)))
    if (length(numeric_index) > 0){
      for (i in 1:length(numeric_index)){
        numeric_index_in_X[i] <- which(attr(X,'assign') == numeric_index[i])
      }
      X[,numeric_index_in_X] <- mean.center(X[,numeric_index_in_X])
      warning("level 1 numeric variables have been mean centered.\n", call. = F, immediate. = T)
    }
    attr(X,'numeric_index') <- numeric_index_in_X
    mf2 <- model.frame(formula = l2_formula, data = data)
    Z_full <- model.matrix(attr(mf2,'terms'), data = mf2)
    if (sum(is.na(Z_full)) > 0) stop('Missing/NA values in level 2 features!')
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
    attr(Z,'varValues') <- get.values(data, mf2, id_index)
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
  result <- list(X = X, Z = Z, Z_full = Z_full, y = y)
  return(result)
}
