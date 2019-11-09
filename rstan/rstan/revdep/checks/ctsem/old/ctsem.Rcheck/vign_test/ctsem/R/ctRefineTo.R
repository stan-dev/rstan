#' ctRefineTo
#' 
#' Fits a ctsem m in a stepwise fashion to help with difficult optimization.
#' 
#' This function fits a sequence of ctsem models increasing in complexity, 
#' starting with a m involving fixed and relatively strong auto effects, no cross effects, no predictors, and no off-diagonal covariances.
#' For many models this can improve the speed and robustness of fitting
#' @param datawide Data in ctsem wide format
#' @param ctmodelobj A continuous time m specified via the \code{\link{ctModel}} function.
#' @param ... additional parameters to pass to \code{\link{ctFit}}.
#' @param modfunc function to run prior to each optimization step, that takes ctsem fit object, modifies it as desired, and returns the fit object.
#' @return Returns a fitted ctsem object in the same manner as \code{\link{ctFit}}.
#' @export

 

ctRefineTo<-function(datawide,ctmodelobj,modfunc=NULL,...){
  
  message('Fitting simplified m with any free DRIFT diagonals fixed to -.3, no free cross effects, no trait / diffusion correlations, no predictors, LAMBDA loadings fixed to 1')

  m<-ctmodelobj
  if(!is.null(m$TRAITVAR) & m$n.latent > 1) {
    m$TRAITVAR[!diag(m$n.latent)]<-0
    m$T0TRAITEFFECT[!diag(m$n.latent)]<-0
  }
  if( m$n.latent > 1) m$DIFFUSION[!diag(m$n.latent)]<-0
  if( m$n.manifest > 1) m$MANIFESTVAR[!diag(m$n.manifest)]<-0
  #if( m$n.latent > 1) m$T0VAR[!diag(m$n.latent)]<-0
  if(!is.null(m$MANIFESTTRAITVAR) & m$n.manifest > 1) m$MANIFESTTRAITVAR[!diag(m$n.manifest)]<-0
  m$DRIFT[row(m$DRIFT) != col(m$DRIFT)][is.na(suppressWarnings(as.numeric(m$DRIFT[row(m$DRIFT) != col(m$DRIFT)])))]<- -.00001
  m$DRIFT[is.na(suppressWarnings(as.numeric(m$DRIFT)))]<-diag(-.3,m$n.latent)[is.na(suppressWarnings(as.numeric(m$DRIFT)))]
  m$LAMBDA[suppressWarnings(is.na(as.numeric(m$LAMBDA)))]<-1
  if(m$n.TDpred > 0) m$TDPREDEFFECT<-matrix(0,nrow=m$n.latent,ncol=m$n.TDpred)
  if(m$n.TDpred > 0) m$T0TDPREDCOV<-matrix(0,nrow=m$n.latent,ncol=(m$n.TDpred*(m$Tpoints)))
  if(m$n.TIpred > 0) m$TIPREDEFFECT<-matrix(0,nrow=m$n.latent,ncol=m$n.TIpred)
  if(m$n.TIpred > 0) m$T0TIPREDEFFECT<-matrix(0,nrow=m$n.latent,ncol=m$n.TIpred)
  if(m$n.TDpred > 0 & m$n.TIpred > 0) m$TDTIPREDCOV<-matrix(0,nrow=(m$n.TDpred*(m$Tpoints)),ncol=m$n.TIpred)
  if(!is.null(m$TRAITTDPREDCOV)) m$TRAITTDPREDCOV<-matrix(0,nrow=m$n.latent,ncol=(m$n.TDpred*(m$Tpoints)))
  if(m$n.TDpred > 0) m$TDPREDVAR[lower.tri(m$TDPREDVAR)]<-0
  if(m$n.TIpred > 0) m$TIPREDVAR[lower.tri(m$TIPREDVAR)]<-0
  fit<-ctFit(datawide,m,nofit=TRUE,...)
  if(!is.null(modfunc)) fit<-modfunc(fit)
  fit$mxobj<-mxTryHard(fit$mxobj, initialTolerance=1e-14,
    initialGradientIterations=1,
    showInits=FALSE, checkHess=TRUE, greenOK=FALSE, 
    iterationSummary=FALSE, bestInitsOutput=FALSE, verbose=fit$ctfitargs$verbose,
    extraTries=fit$ctfitargs$retryattempts, loc=1, scale=0.1, paste=FALSE)
  startValues<-OpenMx::omxGetParameters(fit$mxobj)
  
  
  
    message('Adding correlations, autoeffects, and any predictors, but no free cross effects')
  m<-ctmodelobj
  m$DRIFT[row(m$DRIFT) != col(m$DRIFT)][is.na(suppressWarnings(as.numeric(m$DRIFT[row(m$DRIFT) != col(m$DRIFT)])))]<- -.00001
  fit<-ctFit(datawide,m,nofit=TRUE,omxStartValues=startValues,...)
  if(!is.null(modfunc)) fit<-modfunc(fit)
  fit$mxobj<-mxTryHard(fit$mxobj, initialTolerance=1e-14,
    initialGradientIterations=1,
    showInits=FALSE, checkHess=TRUE, greenOK=FALSE, 
    iterationSummary=FALSE, bestInitsOutput=FALSE, verbose=fit$ctfitargs$verbose,
    extraTries=fit$ctfitargs$retryattempts, loc=1, scale=0.1, paste=FALSE)
    startValues<-OpenMx::omxGetParameters(fit$mxobj)
    oldfit<-fit
  
    if(any(is.na(suppressWarnings(as.numeric(ctmodelobj$DRIFT[!diag(ctmodelobj$n.latent)]))))){
  message('Fitting complete user specified model')
  m<-ctmodelobj
  fit<-ctFit(datawide,m,nofit=TRUE,omxStartValues=startValues,crossEffectNegStarts=FALSE,...)
  if(!is.null(modfunc)) fit<-modfunc(fit)
  fit$mxobj<-mxTryHard(fit$mxobj, initialTolerance=1e-14,
    initialGradientIterations=1,
    showInits=FALSE, checkHess=TRUE, greenOK=FALSE, 
    iterationSummary=FALSE, bestInitsOutput=FALSE, verbose=fit$ctfitargs$verbose,
    extraTries=fit$ctfitargs$retryattempts, loc=1, scale=0.1, paste=FALSE)
    }
  
return(fit)
}
  
  
