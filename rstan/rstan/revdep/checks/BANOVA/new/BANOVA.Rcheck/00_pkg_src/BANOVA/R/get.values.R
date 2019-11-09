get.values <-
function(data, mf, id_index = array(dim = 0), intercept_names = character(0)){
  #
  # Args:
  #   data     : original data frame
  #   mf       : model frame
  #   id_index : corresponding row index for Z matrix values
  # Return:
  #   a list of values included in mf (y values, continuous variable before mean centering-shouldn't be used)
  
# the order of var_names in results matters
  var_names <- c(intercept_names, attr(mf,'names'))
  results <- list()
  if (length(var_names) > 0){
    for (i in 1: length(var_names)){
      if (length(id_index) == 0){
        eval(parse(text = paste('results[[',i,']] <- data$',var_names[i],sep='')))
      }else{
        eval(parse(text = paste('results[[',i,']] <- data$',var_names[i], '[id_index]', sep='')))
      }
      eval(parse(text = paste("attr(results[[",i,"]], 'var_names') <- '", var_names[i],"'",sep="")))
    }
    
  }
  return(results)
}
