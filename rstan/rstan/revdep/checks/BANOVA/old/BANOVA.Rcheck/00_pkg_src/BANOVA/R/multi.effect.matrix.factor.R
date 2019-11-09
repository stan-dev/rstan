multi.effect.matrix.factor <-
function (n_choice, factors, assign = array(dim = 0), index_factor = NA, numeric_index = array(dim = 0), contrast = NULL){

# old version
#  result <- list()
# for (n_c in 1:n_choice){
#   if (length(assign) != 0){
#     index <- which(assign == index_factor)
#     level <- length(index) + 1
#     effect_matrix <- matrix(0, nrow = level, ncol = length(assign))
#     effect_matrix[,index] <- contr.sum(level)
#     if (n_c != 1)
#       effect_matrix[,n_c-1] <- 1 # grand mean included
#     #if (length(numeric_index) > 0) effect_matrix[, numeric_index] <- 1 # consider the covariates effect
#     attr(effect_matrix, 'levels') <- factors
#   }
#   result[[n_c]] <- effect_matrix
# }
  
  # new version
  tmp_contrasts <- getOption("contrasts")
  options(contrasts = rep("contr.sum",2))
  result <- list()
  for (n_c in 1:n_choice){
    if (length(assign) != 0){
      level <- as.factor(levels(factors))
      var_name <- attr(factors,'var_names')
      ### 1.1.2
      level <- assign_contrast_factor(level, var_name, contrast)
      ###
      eval(parse(text = paste(var_name,'<- level', sep = '')))
      # with column names, and include an intercept
      eval(parse(text = paste('effect_matrix <- model.matrix(~',var_name,', data = ', var_name,')', sep='')))
      # Note hard coded here, be consistent with 'multi.design.matrix.R'
      if (n_c > 1){
        intercept_name = paste('choice', n_c,'_Intercept', sep = "")
        temp <- colnames(effect_matrix)
        temp[grep('Intercept', temp)] <- intercept_name
        colnames(effect_matrix) <- temp
      }else{
        # remove choice 1 intercept, identification problem
        #effect_matrix <- effect_matrix[-grep('Intercept', temp)]
      }
      attr(effect_matrix, 'levels') <- levels(factors)
    }else{
      effect_matrix <- NA
    }
    result[[n_c]] <- effect_matrix
  }
  options(contrasts = tmp_contrasts)
  return(result)
}
