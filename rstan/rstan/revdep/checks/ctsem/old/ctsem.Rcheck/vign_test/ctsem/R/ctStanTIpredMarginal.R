#' Plot marginal relationships between covariates and parameters for a ctStanFit object.
#'
#' @param fit ctStanFit object.
#' @param tipred Integer representing which tipred to use -- integer corresponds to TIpredNames specification.
#' @param pars Subject level matrices from the ctStanFit output -- e.g, 'DRIFT' or 'DIFFUSION'.
#' @param probs vector of 3 quantile probabilities, the 2nd will be plotted as a line, the 
#' outer two as shaded regions.
#' @param useimputed Logical, include imputed tipreds or only observed?
#' @param plot Logical, whether to plot.
#'
#' @return If \code{plot=TRUE}, nothing, otherwise an array that can be used with ctPlotArray.
#' @export
#'
#' @examples
#' \dontrun{
#' ctStanTIpredMarginal(ctstantestfit,pars='CINT',tipred=3)
#' }
ctStanTIpredMarginal<-function(fit,tipred,pars,probs=c(.025,.5,.975),useimputed=TRUE, plot=TRUE){
  e<-extract(fit)
  p<- e[[pars]]
  if(useimputed) tipreds <- ctCollapse(e$tipreds,1,mean) else tipreds <- fit$data$tipredsdata
  p<-plyr::aaply(p,2:length(dim(p)),quantile,probs=probs,.drop=FALSE)
  pdim<-dim(p)
  p<-array(p,c(pdim[1],prod(pdim[c(-1,-length(pdim))]),pdim[length(pdim)]))
  for(dimi in 2:(length(pdim)-1)){
    if(dimi==2) pnames <- 1:pdim[dimi] else pnames <- paste0(pnames,', ',rep(1:pdim[dimi],each=length(pnames)))
  }
  dimnames(p) <- list(fit$setup$idmap[,'original'],
    paste0(pars,'[',pnames,']'),
    paste0('Q',probs*100,'%'))
  
  #sort
  p <- p[order(tipreds[,tipred[1],drop=FALSE]), ,,drop=FALSE]
  tipreds <- tipreds[order(tipreds[,tipred[1],drop=FALSE]),,drop=FALSE]
  colnames(tipreds) <- fit$ctstanmodel$TIpredNames
  input <- list(y=p,x=tipreds[,tipred[1],drop=FALSE])
  
  if(plot){
    for(ci in 1:(dim(p)[2])){
      intemp<-input
      intemp$y <- intemp$y[,ci, ,drop=FALSE]
      ctPlotArray(intemp)
      # plot(rep(tipreds[,tipred],pdim[1]),p[,ci,drop=FALSE],xlab=fit$ctstanmodel$TIpredNames[tipred],ylab=colnames(p)[ci],pch=pch,pcol=pcol)
    }
  }
}
