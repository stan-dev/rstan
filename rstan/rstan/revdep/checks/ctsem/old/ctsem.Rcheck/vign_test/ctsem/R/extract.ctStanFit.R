#' Extract samples from a ctStanFit object
#'
#' @param object ctStanFit object, samples may be from Stan's HMC, or the importance sampling approach of ctsem.
#' @param ... additional arguments to pass to \code{rstan::extract}.
#' @return Array of posterior samples.
#' @aliases extract
#' @examples
#' e = extract(ctstantestfit)
#' head(e)
#' @export
extract <- function(object,...){
  if(!class(object) %in% c('ctStanFit', 'stanfit')) stop('Not a ctStanFit or stanfit object')
  if(class(object)=='stanfit') out <- rstan::extract(object,...) else{
  if(class(object$stanfit)=='stanfit') out <- rstan::extract(object$stanfit,...)
  if(class(object$stanfit)!='stanfit') out <- object$stanfit$transformedpars
  out$Ygen[out$Ygen==99999] <- NA
  }
  return(out)
}
