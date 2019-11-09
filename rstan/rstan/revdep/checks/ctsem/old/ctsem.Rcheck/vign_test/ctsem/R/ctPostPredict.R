#' Posterior predictive type check for ctsemFit.
#' 
#' Samples data according to the ctsemFit object, computes quantiles over time based on model fit,
#' plots these against original data.
#'
#' @param fit object of class ctsemFit as returned from \code{\link{ctFit}}
#' @param timestep positive value denoting the time interval to use for sampling.
#' @param n.subjects Number of subjects worth of data to sample.
#' @param probs Vector of values between 0 and 1 denoting quantiles to generate. 
#' For plotting, vector should be of length 3 and values should be rising.
#' @param plot Whether to plot or return the generated data.
#' @param ctPlotArrayArgs additional arguments to pass to \code{\link{ctPlotArray}} function,
#' for plotting generated distributions.
#' @param indPlotArgs list of parameters to pass to ctIndplot, for 
#' plotting original data. Only used if plot=TRUE.
#' @param mfrow 2 dimensional integer vector defining number of rows and columns of plots,
#' as per the mfrow argument to \code{\link[graphics]{par}}.
#' 'auto' determines automatically, to a maximum of 4 by 4, while \code{NULL} 
#' uses the current system setting.
#'
#' @return Either nothing (if plot=TRUE) or an array containing generated data over quantiles.
#' @export
#'
#' @examples
#' data("AnomAuth")
#' AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
#'   Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = 'auto') 
#' AnomAuthFit <- ctFit(AnomAuth, AnomAuthmodel)
#' ctPostPredict(AnomAuthFit,timestep=.5,n.subjects=100)

ctPostPredict<-function(fit,timestep=.1,n.subjects = 100,probs=c(.025,.5,.975),
  plot=TRUE, ctPlotArrayArgs=list(grid=FALSE,legend=FALSE), 
  indPlotArgs=list(colourby = 'subject',lwd=2,new=FALSE,type='p',opacity=.3),
  mfrow='auto'){

pcheck=ctGenerateFromFit(fit = fit,n.subjects = n.subjects,wide=FALSE,timestep = timestep)
timeupper=max(pcheck[,'time'])
pcheck=array(t(pcheck),dim=c(ncol(pcheck),nrow(pcheck)/n.subjects,n.subjects),
  dimnames = list(vars=colnames(pcheck),
    time=paste0('T',0:(nrow(pcheck)/n.subjects-1)),
    id=paste0('id_',1:n.subjects)))
pcheck<-plyr::aaply(probs,1,function(x) 
  ctCollapse(pcheck[fit$ctmodelobj$manifestNames,,,drop=FALSE],3,quantile,probs=x),.drop=FALSE)
dimnames(pcheck)=list(quantile=probs,vars=dimnames(pcheck)[[2]],time=dimnames(pcheck)[[3]])
pcheck=aperm(pcheck,c(3,2,1))




#plotting

if(plot==TRUE){
  
  paroriginal<-graphics::par()[c('mfrow')]
  
  if(!is.null(mfrow)){
    if(mfrow=='auto') {
      graphics::par(mfrow=c(min(3,grDevices::n2mfrow( fit$ctmodelobj$n.manifest)[2]), 
        min(3,n2mfrow( fit$ctmodelobj$n.manifest)[1])))
    }
    if(mfrow!='auto') graphics::par(mfrow=mfrow)
  }
  
  #convert kalman data to wide
  if(is.null(fit$mxobj$expectation$P0))  data<-fit$mxobj$data$observed #if not kalman based fit
  if(!is.null(fit$mxobj$expectation$P0)) { 
    data<- suppressMessages(ctLongToWide(fit$mxobj$data$observed,id='id',time='dT1',
      manifestNames=fit$ctmodelobj$manifestNames,
      TDpredNames=fit$ctmodelobj$TDpredNames,
      TIpredNames=fit$ctmodelobj$TIpredNames))
    data<-data[,-which(colnames(data)=='T0'),drop=FALSE]
    colnames(data)[which(colnames(data) %in% paste0('T',1:(fit$ctmodelobj$Tpoints-1)))]<-paste0('dT',1:(fit$ctmodelobj$Tpoints-1))
  }
  
  
  
 for(i in 1:fit$ctmodelobj$n.manifest){
  ctPlotArrayArgs$input$y = pcheck[,i,,drop=FALSE]
  ctPlotArrayArgs$input$x = seq(0,timeupper,timestep)
  ctPlotArrayArgs$plotcontrol=list(ylab='Values',xlab='Time', main=fit$ctmodelobj$manifestNames[i],
    ylim=range(c(ctPlotArrayArgs$yarray,data[,paste0(fit$ctmodelobj$manifestNames[i],'_T',0:(fit$ctmodelobj$Tpoints-1))]),na.rm=TRUE))
  do.call(ctPlotArray,ctPlotArrayArgs)

  
  indPlotArgs$datawide = data
  indPlotArgs$n.subjects = nrow(data)
  indPlotArgs$vars=i
  indPlotArgs$n.manifest = fit$ctmodelobj$n.manifest
  indPlotArgs$Tpoints = fit$ctmodelobj$Tpoints
  indPlotArgs$legend=FALSE
  
suppressMessages(do.call(ctIndplot,indPlotArgs))
}
do.call(graphics::par,paroriginal) #end plotting and reset graphics
} 

if(plot!=TRUE) return(pcheck)
}
