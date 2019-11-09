effect.matrix.interaction <-
function (interaction_factors = list(), 
          assign = array(dim = 0), 
          main_eff_index = array(dim = 0), 
          index_inter_factor = NA, 
          numeric_index = array(dim = 0),
          contrast = NULL){
# generate the effect matrix for each interaction, numerical covariates excluded, TODO: may have various levels for numeric variables
# Args:
#     interaction_factors       : a list of values of factors in the interaction to generate the effect matrix
#     assign        : Not used for the new method, to remove
#     index_factor  : Not used, to remove
#     numeric_index : Not used, to remove
# Returns:
#     an effect matrix
#
# old method, hard coded
if(0){
  # two way interactions
  if (length(interaction_factors) == 2){
    factor1 <- interaction_factors[[1]]
    factor2 <- interaction_factors[[2]]
    n1 <- length(factor1)
    n2 <- length(factor2)
    temp1 <- array(NA,dim = c(n1*n2,1))
    temp2 <- array(NA,dim = c(n1*n2,1))
    for (i in 1:n1)
      for (j in 1: n2){
        temp1[(i-1)*n2 + j] <- factor1[i]
        temp2[(i-1)*n2 + j] <- factor2[j]
      }
    tempdata <- data.frame(cbind(temp1,temp2)) # everything converted to 1,2,3...
    tempdata[,1] <- as.factor(tempdata[,1])
    tempdata[,2] <- as.factor(tempdata[,2])
    tempmatrix <- model.matrix(~ temp1 * temp2, data = tempdata)
    tempIntermatrix <- tempmatrix[,-(1:(n1-1+n2-1+1))] # intercept:first column is excluded
    tempMainmatrix <- tempmatrix[,2:(n1-1+n2-1+1)]
    index <- which(assign == index_inter_factor)
    effect_matrix <- matrix(0, nrow = n1*n2, ncol = length(assign))
    effect_matrix[,index] <- tempIntermatrix
    effect_matrix[,1] <- 1 # grand mean included
    #if (length(numeric_index) > 0) effect_matrix[, numeric_index] <- 1 # consider the covariates effect
    #include main effects
    effect_matrix[, c(which(assign == main_eff_index[1]), which(assign == main_eff_index[2]))] <- tempMainmatrix
    attr(effect_matrix, 'levels') <- cbind(temp1,temp2)
  }else{
    # multi way interactions
    effect_matrix <- NA
  }
  return(effect_matrix)
}

# new method including higher number (>= 2) of interactions
  tmp_contrasts <- getOption("contrasts")
  options(contrasts = rep("contr.sum",2))
  if (length(interaction_factors) >= 2){
    num_inter = length(interaction_factors)
    levels_inter = list()
    names_inter = list()
    for (i in 1:num_inter){
      names_inter[[i]] = attr(interaction_factors[[i]], 'var_names')
      levels_inter[[attr(interaction_factors[[i]], 'var_names')]] = levels(interaction_factors[[i]])
    }
    factors_inter = expand.grid(levels_inter)
    factors_inter_factor = as.data.frame(lapply(factors_inter, as.factor))
    #### 1.1.2
    factors_inter_factor <- assign_contrast(factors_inter_factor, contrast)
    ####
    formula_inter = paste(names_inter, collapse = '*')
    eval(parse(text = paste('effect_matrix <- model.matrix(~',formula_inter,', data = factors_inter_factor)', sep='')))
    attr(effect_matrix, 'levels') = as.matrix(factors_inter)
  }else{
    effect_matrix = NA
  }
  options(contrasts = tmp_contrasts)
  return(effect_matrix)
}
