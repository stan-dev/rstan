#' Extract functions of multiple variables from a stanfit object
#'
#' Can be useful for determining quantiles or plotting multidimensional regions -- 
#' for instance in case of colinearity of predictors.
#' @param stanfit object of class stanfit.
#' @param parstrings vector of strings containing partial (or full) matches of parameter names.
#' When more than one string is passed, functions are computed based on the combination of the first match for
#' each string, then the second match for each string, etc. The first match of the first string is only ever combined
#' with the first match of the second, similarly for the 2nd match, etc.
#' @param prefuncstring string containing front element of function. E.g., 'exp(' for an exponential region.
#' @param joinfuncstring string used to join the (possibly) multiple parameters involved.
#' @param postfuncstring string containing end element of function. E.g., ') *2' to multiply the result by 2.
#'
#' @return matrix of values of the specified interactions at each iteration. 
#' @export
#'
#' @examples
#' temp<-stan_confidenceRegion(stanfit=ctstantestfit$stanfit, 
#'   parstrings=c('pop_DRIFT[1,2]','pop_DRIFT[2,1]'))
#' t(apply(temp,2,quantile))
stan_confidenceRegion <-function(stanfit,parstrings,prefuncstring='(', joinfuncstring=' + ',postfuncstring=')'){
  mc=As.mcmc.list(stanfit)
  mc=do.call(rbind,mc)
  
  pars <- lapply(parstrings,function(x) paste0(colnames(mc)[grep(x,colnames(mc),fixed=TRUE)]))
  parsref <- lapply(parstrings,function(x) paste0('mc[,"',colnames(mc)[grep(x,colnames(mc),fixed=TRUE)],'"]'))
  
  # if(length(parstrings) > 1 & matchingindices==TRUE & !all(lapply(pars,length)==length(pars[[1]]))) stop ('matchingindices=TRUE but unequal numbers of parameters found matching parstrings')
  
  for(pari in 1:length(pars[[1]])){
    a <- cbind(eval(parse(text=paste0(prefuncstring, paste0(lapply(parsref,function(x) x[pari]),collapse=joinfuncstring),postfuncstring))))
    colnames(a) <- paste0(prefuncstring, paste0(lapply(pars,function(x) x[pari]),collapse=joinfuncstring),postfuncstring)
    if(pari==1) out<-a else out <- cbind(out,a)
  }
  
  return(out)
}
