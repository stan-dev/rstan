#' Compute functions of matrices from samples of a stanfit object
#'
#' @param stanfit object of class stanfit.
#' @param object name of stan sub object from stanfit to use for calculations.
#' @param objectindices matrix of indices, with the number of columns matching 
#' the number of dimensions of the object. 'all' computes \code{which( array(1,objdims)==1,arr.ind=TRUE)},
#' where objdims is what would be returned by dim(object) if the object existed in the R environment.
#' @param calc string containing R calculation to evaluate, with the string 'object' in place of the actual object name.
#' @param summary if FALSE, a iterations * parameters matrix is returned, if TRUE, 
#' rstan::monitor is first run on the output.
#'
#' @return matrix of values of the specified interactions at each iteration. 
#' @export
#'
#' @examples
#' temp<-stan_postcalc(stanfit=ctstantestfit$stanfit, 
#'   object='DRIFT', objectindices='all', calc='exp(object)')
stan_postcalc <-function(stanfit,object,calc='object', objectindices='all', summary=TRUE){
  mc=As.mcmc.list(stanfit)
  m=do.call(rbind,mc)
  outdims = dim(stanfit@inits[[1]][[object]]) #complete
  if(objectindices=='all') objectindices <- which( array(1,outdims)==1,arr.ind=TRUE) else {
  if(ncol(as.matrix(objectindices))!= length(outdims)) stop('Number of columns of object indices must match number of dimensions in object')
  }
  objectindices = objectindices[,outdims!=1] #subset
  outdims = outdims[outdims!=1] #subset

  out=array(apply(m,1,function(x) {
    object = array(array(relist(x,skeleton = stanfit@inits[[1]])[[object]],dim=outdims)[objectindices],outdims)
    ret = eval(parse(text=calc))
  } ),dim = c(outdims,nrow(m)) )
  
  if(summary) {
    parnames <- array(1,dim=outdims)
    parnames <- which(parnames==1,arr.ind = TRUE)
    parnames <- apply(parnames,1,function(x) paste0(x, collapse=', '))
    parnames <- paste0(object,'[',parnames,']')
    out <- array(out,dim=c(prod(outdims),  nrow(mc[[1]]), length(mc)))

    out <- aperm(out,c(2,3,1))
    dimnames(out)[[3]] <- parnames
    out <- monitor(out,warmup=0,digits_summary = 2)
  }
  
  return(out)
}

