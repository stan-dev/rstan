BANOVA.TNormal <-
function(l1_formula = 'NA', l2_formula = 'NA', data, id, l1_hyper, l2_hyper, burnin, sample, thin, adapt, conv_speedup, jags){
  cat('Model initializing...\n')
  if (l1_formula == 'NA'){
    stop("Formula in level 1 is missing or not correct!")
  }else{
    mf1 <- model.frame(formula = l1_formula, data = data)
    # check y, if it is numeric
    y <- model.response(mf1)
    if (class(y) != 'numeric'){
      warning("The response variable must be numeric (data class also must be 'numeric')..")
      y <- as.numeric(y)
      warning("The response variable has been converted to numeric..")
    }
  }
  single_level = F
  if (l2_formula == 'NA'){
    # one level models
    single_level = T
    # check each column in the dataframe should have the class 'factor' or 'numeric', no other classes such as 'matrix'...
    for (i in 1:ncol(data)){
      if(class(data[,i]) != 'factor' && class(data[,i]) != 'numeric' && class(data[,i]) != 'integer') stop("data class must be 'factor', 'numeric' or 'integer'")
      # checking numerical predictors, converted to categorical variables if the number of levels is <= 3
      if ((class(data[,i]) == 'numeric' | class(data[,i]) == 'integer') & length(unique(data[,i])) <= 3){
        data[,i] <- as.factor(data[,i])
        warning("Variables(levels <= 3) have been converted to factors")
      }
    }
    n <- nrow(data)
    uni_id <- unique(id)
    num_id <- length(uni_id)
    new_id <- rep(0, length(id)) # store the new id from 1,2,3,...
    for (i in 1:length(id))
      new_id[i] <- which(uni_id == id[i])
    id <- new_id
    dMatrice <- design.matrix(l1_formula, l2_formula, data = data, id = id)
    JAGS.model <- JAGSgen.TNormal(dMatrice$X, dMatrice$Z, l1_hyper, l2_hyper, conv_speedup)
    JAGS.data <- dump.format(list(n = n, y = dMatrice$y, X = dMatrice$X))
    result <- run.jags (model = JAGS.model$sModel, data = JAGS.data, inits = JAGS.model$inits, n.chains = 1,
                        monitor = c(JAGS.model$monitorl1.parameters, JAGS.model$monitorl2.parameters), 
                        burnin = burnin, sample = sample, thin = thin, adapt = adapt, jags = jags, summarise = FALSE, 
                        method="rjags")
    samples <- result$mcmc[[1]]
    # find the correct samples, in case the order of monitors is shuffled by JAGS.
    n_p_l1 <- length(JAGS.model$monitorl1.parameters)
    index_l1_param<- array(0,dim = c(n_p_l1,1))
    for (i in 1:n_p_l1)
      index_l1_param[i] <- which(colnames(result$mcmc[[1]]) == JAGS.model$monitorl1.parameters[i])
    if (length(index_l1_param) > 1)
      samples_l1_param <- result$mcmc[[1]][,index_l1_param]
    else
      samples_l1_param <- matrix(result$mcmc[[1]][,index_l1_param], ncol = 1)
    colnames(samples_l1_param) <- colnames(result$mcmc[[1]])[index_l1_param]
    cat('Constructing ANOVA/ANCOVA tables...\n')
    dMatrice$Z <-  array(1, dim = c(1,1), dimnames = list(NULL, ' '))
    attr(dMatrice$Z, 'assign') <- 0
    attr(dMatrice$Z, 'varNames') <- " "
    samples_l2_param <- NULL
    anova.table <- table.ANCOVA(samples_l2_param, dMatrice$Z, dMatrice$X, samples_l1_param, array(y, dim = c(length(y), 1))) # for ancova models
    coef.tables <- table.coefficients(samples_l1_param, JAGS.model$monitorl1.parameters, colnames(dMatrice$Z), colnames(dMatrice$X), 
                                      attr(dMatrice$Z, 'assign') + 1, attr(dMatrice$X, 'assign') + 1)
    pvalue.table <- table.pvalue(coef.tables$coeff_table, coef.tables$row_indices, l1_names = attr(dMatrice$Z, 'varNames'), 
                                 l2_names = attr(dMatrice$X, 'varNames'))
    conv <- conv.geweke.heidel(samples_l1_param, colnames(dMatrice$Z), colnames(dMatrice$X))
    mf2 <- NULL
    class(conv) <- 'conv.diag'
    cat('Done...\n')
    
  }else{
    mf2 <- model.frame(formula = l2_formula, data = data)
    # check each column in the dataframe should have the class 'factor' or 'numeric', no other classes such as 'matrix'...
    for (i in 1:ncol(data)){
      if(class(data[,i]) != 'factor' && class(data[,i]) != 'numeric' && class(data[,i]) != 'integer') stop("data class must be 'factor', 'numeric' or 'integer'")
      # response_name <- attr(mf1,"names")[attr(attr(mf1, "terms"),"response")]
      # checking missing predictors, already checked in design matrix
      # if(i != which(colnames(data) == response_name) & sum(is.na(data[,i])) > 0) stop("Data type error, NAs/missing values included in independent variables") 
      #if(i != which(colnames(data) == response_name) & class(data[,i]) == 'numeric')
      #  data[,i] = data[,i] - mean(data[,i])
      # checking numerical predictors, converted to categorical variables if the number of levels is <= 3
      if ((class(data[,i]) == 'numeric' | class(data[,i]) == 'integer') & length(unique(data[,i])) <= 3){
        data[,i] <- as.factor(data[,i])
        warning("Variables(levels <= 3) have been converted to factors")
      }
    }
    n <- nrow(data)
    uni_id <- unique(id)
    num_id <- length(uni_id)
    new_id <- rep(0, length(id)) # store the new id from 1,2,3,...
    for (i in 1:length(id))
      new_id[i] <- which(uni_id == id[i])
    id <- new_id
    dMatrice <- design.matrix(l1_formula, l2_formula, data = data, id = id)
    JAGS.model <- JAGSgen.TNormal(dMatrice$X, dMatrice$Z, l1_hyper, l2_hyper, conv_speedup)
    JAGS.data <- dump.format(list(n = n, id = id, M = num_id, y = dMatrice$y, X = dMatrice$X, Z = dMatrice$Z))
    result <- run.jags (model = JAGS.model$sModel, data = JAGS.data, inits = JAGS.model$inits, n.chains = 1,
                        monitor = c(JAGS.model$monitorl1.parameters, JAGS.model$monitorl2.parameters), 
                        burnin = burnin, sample = sample, thin = thin, adapt = adapt, jags = jags, summarise = FALSE, 
                        method="rjags")
    samples <- result$mcmc[[1]]
    # find the correct samples, in case the order of monitors is shuffled by JAGS
    n_p_l2 <- length(JAGS.model$monitorl2.parameters)
    index_l2_param<- array(0,dim = c(n_p_l2,1))
    for (i in 1:n_p_l2)
      index_l2_param[i] <- which(colnames(result$mcmc[[1]]) == JAGS.model$monitorl2.parameters[i])
    if (length(index_l2_param) > 1)
      samples_l2_param <- result$mcmc[[1]][,index_l2_param]
    else
      samples_l2_param <- matrix(result$mcmc[[1]][,index_l2_param], ncol = 1)
    colnames(samples_l2_param) <- colnames(result$mcmc[[1]])[index_l2_param]
    
    n_p_l1 <- length(JAGS.model$monitorl1.parameters)
    index_l1_param<- array(0,dim = c(n_p_l1,1))
    for (i in 1:n_p_l1)
      index_l1_param[i] <- which(colnames(result$mcmc[[1]]) == JAGS.model$monitorl1.parameters[i])
    if (length(index_l1_param) > 1)
      samples_l1_param <- result$mcmc[[1]][,index_l1_param]
    else
      samples_l1_param <- matrix(result$mcmc[[1]][,index_l1_param], ncol = 1)
    colnames(samples_l1_param) <- colnames(result$mcmc[[1]])[index_l1_param]
    
    #anova.table <- table.ANOVA(samples_l1_param, dMatrice$X, dMatrice$Z)
    cat('Constructing ANOVA/ANCOVA tables...\n')
    anova.table <- table.ANCOVA(samples_l1_param, dMatrice$X, dMatrice$Z, samples_l2_param) # for ancova models
    coef.tables <- table.coefficients(samples_l2_param, JAGS.model$monitorl2.parameters, colnames(dMatrice$X), colnames(dMatrice$Z), 
                                      attr(dMatrice$X, 'assign') + 1, attr(dMatrice$Z, 'assign') + 1)
    pvalue.table <- table.pvalue(coef.tables$coeff_table, coef.tables$row_indices, l1_names = attr(dMatrice$X, 'varNames'), 
                                 l2_names = attr(dMatrice$Z, 'varNames'))
    conv <- conv.geweke.heidel(samples_l2_param, colnames(dMatrice$X), colnames(dMatrice$Z))
    class(conv) <- 'conv.diag'
    cat('Done...\n')
  }
  return(list(anova.table = anova.table,
              coef.tables = coef.tables,
              pvalue.table = pvalue.table, 
              conv = conv,
              dMatrice = dMatrice, samples_l1_param = samples_l1_param, 
              samples_l2_param = samples_l2_param, data = data, mf1 = mf1, mf2 = mf2,
              JAGSmodel = JAGS.model$sModel, single_level = single_level, model_name = "BANOVA.T"))
}
