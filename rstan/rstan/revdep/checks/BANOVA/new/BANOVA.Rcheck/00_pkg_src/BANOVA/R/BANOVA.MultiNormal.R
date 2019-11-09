BANOVA.MultiNormal <-
function(l1_formula = 'NA', l2_formula = 'NA', dataX, dataZ, y, id, l2_hyper, burnin, sample, thin, adapt, conv_speedup, jags){
  cat('Model initializing...\n')
  # TODO one level model
  if (l2_formula == 'NA')
    stop("The level 2 formula is not specified, please check BANOVA.run for single level models.")
  # check y, if it is integers
  if (class(y) != 'integer'){
    warning("The response variable must be integers (data class also must be 'integer')..")
    y <- as.integer(as.character(y))
    warning("The response variable has been converted to integers..")
  }
  DV_sort <- sort(unique(y))
  n_categories <- length(DV_sort)
  if (n_categories < 2) stop('The number of categories must be greater than 1!')
  if (DV_sort[1] != 1 || DV_sort[n_categories] != n_categories) stop('Check if response variable follows categorical distribution!') 
  # check each column in the dataframe should have the class 'factor' or 'numeric', no other classes such as 'matrix'...
  for (i in 1:ncol(dataZ)){
    if(class(dataZ[,i]) != 'factor' && class(dataZ[,i]) != 'numeric' && class(dataZ[,i]) != 'integer') stop("data class must be 'factor', 'numeric' or 'integer'")
    # checking missing predictors, already checked in design matrix
    # if(sum(is.na(dataZ[,i])) > 0) stop("Data type error, NAs/missing values included in independent variables") 
    #if(class(dataZ[,i]) == 'numeric')
    #  dataZ[,i] = dataZ[,i] - mean(dataZ[,i])
    # checking numerical predictors, converted to categorical variables if the number of levels is <= 3
    if ((class(dataZ[,i]) == 'numeric' | class(dataZ[,i]) == 'integer') & length(unique(dataZ[,i])) <= 3){
      dataZ[,i] <- as.factor(dataZ[,i])
      warning("Between-subject variables(levels <= 3) have been converted to factors")
    }
  }
  for (i in 1:length(dataX))
    for (j in 1:ncol(dataX[[i]])){
      if(class(dataX[[i]][,j]) != 'factor' && class(dataX[[i]][,j]) != 'numeric' && class(dataX[[i]][,j]) != 'integer') stop("data class must be 'factor', 'numeric' or 'integer'")
      # checking missing predictors, already checked in design matrix
      # if(sum(is.na(dataX[[i]][,j])) > 0) stop("Data type error, NAs/missing values included in independent variables") 
      #if(class(dataX[[i]][,j]) == 'numeric')
      #  dataX[[i]][,j] = dataX[[i]][,j] - mean(dataX[[i]][,j])
      # checking numerical predictors, converted to categorical variables if the number of levels is <= 3
      if ((class(dataX[[i]][,j]) == 'numeric' | class(dataX[[i]][,j]) == 'integer') & length(unique(dataX[[i]][,j])) <= 3){
        dataX[[i]][,j] <- as.factor(dataX[[i]][,j])
        warning("Within-subject variables(levels <= 3) have been converted to factors")
      }
    }
  n <- nrow(dataZ)
  uni_id <- unique(id)
  num_id <- length(uni_id)
  new_id <- rep(0, length(id)) # store the new id from 1,2,3,...
  for (i in 1:length(id))
    new_id[i] <- which(uni_id == id[i])
  id <- new_id
  dMatrice <- multi.design.matrix(l1_formula, l2_formula, dataX = dataX, dataZ = dataZ, id = id)
  
  # create 3-dimensional matrix of X for BUGS data
  X_new <- array(0, dim = c(n, ncol(dMatrice$X_full[[1]]), n_categories))
  for (i in 1:n_categories){
    X_new[,,i] <- dMatrice$X_full[[i]]
  }
  JAGS.model <- JAGSgen.multiNormal(X_new, dMatrice$Z, l2_hyper, conv_speedup)
  JAGS.data <- dump.format(list(n = n, id = id, M = num_id, y = y, X = X_new, Z = dMatrice$Z, n_choice = n_categories))
  result <- run.jags (model = JAGS.model$sModel, data = JAGS.data, inits = JAGS.model$inits, n.chains = 1,
                      monitor = c(JAGS.model$monitorl1.parameters, JAGS.model$monitorl2.parameters, JAGS.model$monitor.cutp), 
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
  
  #anova.table <- table.ANOVA(samples_l1_param, dMatrice$X_full[[1]], dMatrice$Z) # only need the colnames of X
  cat('Constructing ANOVA/ANCOVA tables...\n')
  anova.table <- table.ANCOVA(samples_l1_param, dMatrice$X_full[[1]], dMatrice$Z, samples_l2_param) # for ancova models
  coef.tables <- table.coefficients(samples_l2_param, JAGS.model$monitorl2.parameters, colnames(dMatrice$X_full[[1]]), colnames(dMatrice$Z), 
                                    attr(dMatrice$X_full[[1]], 'assign'), attr(dMatrice$Z, 'assign') + 1)
  pvalue.table <- table.pvalue(coef.tables$coeff_table, coef.tables$row_indices, l1_names = attr(dMatrice$X_full[[1]], 'varNames'), 
                               l2_names = attr(dMatrice$Z, 'varNames'))
  conv <- conv.geweke.heidel(samples_l2_param, colnames(dMatrice$X_full[[1]]), colnames(dMatrice$Z))
  class(conv) <- 'conv.diag'
  mf1 <- model.frame(formula = l1_formula, data = dataX[[1]])
  mf2 <- model.frame(formula = l2_formula, data = dataZ)
  cat('Done...\n')
  return(list(anova.table = anova.table,
              coef.tables = coef.tables,
              pvalue.table = pvalue.table, 
              conv = conv,
              dMatrice = dMatrice, samples_l2_param = samples_l2_param, dataX = dataX, dataZ = dataZ,
              mf1 = mf1, mf2 = mf2, n_categories = n_categories,JAGSmodel = JAGS.model$sModel, single_level = F, model_name = "BANOVA.Multinomial"))
}
