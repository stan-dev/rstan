#' Update an already compiled and fit ctStanFit object
#' 
#' Allows one to change data and or model elements that don't require recompiling, then re fit.
#'
#' @param fit ctStanFit object
#' @param datalong data as normally passed to \code{\link{ctStanFit}}
#' @param ctstanmodel model as normally passed to \code{\link{ctStanFit}}
#' @param ... extra args for \code{\link{ctStanFit}}
#' @examples
#' \dontrun{
#'  newm<-ctModel(type='stanct',
#'   n.latent=ctstantestfit$ctstanmodel$n.latent,
#'   n.TDpred=ctstantestfit$ctstanmodel$n.TDpred,
#'   n.TIpred=ctstantestfit$ctstanmodel$n.TIpred,
#'   MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
#'   MANIFESTMEANS=matrix(0,nrow=ctstantestfit$ctstanmodel$n.manifest),
#'   CINT=matrix(c(0,'cint2'),ncol=1),
#'   n.manifest=ctstantestfit$ctstanmodel$n.manifest,
#'   LAMBDA=diag(2))
#'   
#'  newdat <- ctstantestdat
#'  newdat <- newdat[newdat[,'id']!=1,]
#'  newfit <- ctStanUpdModel(ctstantestfit, newdat, newm)
#'  }

ctStanUpdModel <- function(fit, datalong, ctstanmodel,...){
  
  new <-ctStanFit(datalong = datalong, ctstanmodel = ctstanmodel,fit=FALSE,...)
  
  fit$standata <- new$standata
  fit$data <- new$data
  fit$setup <- new$setup
  fit$args <- match.call
  return(fit)
}
