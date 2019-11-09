#' Convert samples from a stanfit object to the unconstrained scale
#'
#' @param fit stanfit object.
#' @param standata only necessary if R session has been restarted since fitting model -- used to reinitialize 
#' stanfit object.
#'
#' @return Matrix containing columns of unconstrained parameters for each post-warmup iteration.
#' @export
#'
#' @examples
#' \dontrun{
#' umat <- stan_unconstrainsamples(ctstantestfit$stanfit, ctstantestfit$standata)
#' }
stan_unconstrainsamples <- function(fit, standata=NA){
  if(class(fit)!='stanfit') stop('not a stanfit object')
  npars <- try(get_num_upars(fit),silent=TRUE) #$stanmodel)
  
  if(class(npars)=='try-error'){ #in case R has been restarted or similar
    if(any(!is.na(standata))){
      newfit <- stan_reinitsf(fit@stanmodel,standata) 
    } 
    else stop('stanfit object must be reinitialized but no data is provided')
  } else newfit <- fit #no need for reinit
  
  cmat=as.matrix(fit)
  # clist=apply(cmat,1,function(x) relist(flesh = x,skeleton = fit@inits[[1]]))
  
  clist=apply(cmat,1,function(x) relistarrays(flesh=x,skeleton=fit@inits[[1]]))

  ulist=matrix(unlist(lapply(clist,function(x) unconstrain_pars(newfit,x))),ncol=length(clist))
  return(ulist)
}
