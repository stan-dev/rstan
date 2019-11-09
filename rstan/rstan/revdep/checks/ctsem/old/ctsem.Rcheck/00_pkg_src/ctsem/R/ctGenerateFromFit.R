#' Generates data according to the model estimated in a ctsemFit object.#' 
#'
#' @param fit object of class ctsemFit as returned from \code{\link{ctFit}}
#' @param timestep positive numeric value indicating the time interval to use for
#' data generation.
#' @param timerange either 'asdata' to calculate range based on data in fit object,
#' or vector of length 2 specifying min and max times for generation.
#' @param n.subjects integer. Number of subjects worth of data to generate
#' @param predictorSubjects vector of integers, or string 'all', defining which 
#' subjects to sample time dependent and independent predictors from.
#' @param ... parameters to pass to ctGenerate function, such as wide=FALSE.
#'
#' @return matrix of generated data
#' @export
#'
#' @examples
#' 
#' data(AnomAuth) 
#' AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
#'   Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2)) 
#' AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)
#' 
#' dwide <- ctGenerateFromFit(AnomAuthfit,timestep=1,n.subjects=5)
#' 
#' par(mfrow=c(1,2))
#' ctIndplot(datawide = dwide,n.subjects = 5,n.manifest = 2,vars=1,Tpoints = 4)
#' ctIndplot(datawide = AnomAuth+rnorm(length(AnomAuth)),vars=1,n.subjects = 5,
#' n.manifest = 2,Tpoints = 4)
#' 
ctGenerateFromFit<-function(fit,timestep='asdata',n.subjects=100,timerange='asdata',
  predictorSubjects='all',...){


gm=ctModelFromFit(fit)

if(!is.null(fit$mxobj$expectation$P0)) { #if fit with kalman filter then data needs rearranging
  dat=suppressMessages(ctLongToWide(datalong = fit$mxobj$data$observed,
  id = 'id',time = 'dT1',manifestNames = gm$manifestNames,TDpredNames = gm$TDpredNames,
  TIpredNames = gm$TIpredNames))
  
  dat=dat[,-which(colnames(dat)=='T0')]
  colnames(dat)[colnames(dat) %in% c(paste0('T',1:(gm$Tpoints-1)))]=paste0('dT',1:(gm$Tpoints-1))
} else dat=fit$mxobj$data$observed

if(predictorSubjects=='all') predictorSubjects=1:(nrow(dat))

if(timerange=='asdata') timerange=c(0,max(apply(dat[,paste0('dT',1:(gm$Tpoints-1))],1,sum,na.rm=TRUE)))

if(timestep!='asdata'){
  if(is.na(as.numeric(timestep)) || as.numeric(timestep) <= 0) stop('timestep must be a positive value or "asdata"')
  gm$Tpoints=length(seq(timerange[1],timerange[2],timestep))
}


out=c()
for(i in 1:n.subjects){
  
  if(gm$n.TDpred + gm$n.TIpred > 0){ #then discretise data for this subject to timestep so that predictor information is accurate
    predSub=sample(x = predictorSubjects,size = 1)
    
    ndlong <- suppressMessages(ctWideToLong(datawide=dat[predSub,,drop=FALSE],Tpoints=fit$ctmodelobj$Tpoints,
      n.manifest=gm$n.manifest,n.TDpred=gm$n.TDpred,n.TIpred=gm$n.TIpred,
      manifestNames=gm$manifestNames,TDpredNames=gm$TDpredNames,TIpredNames=gm$TIpredNames))
    
    ndlong <- suppressMessages(ctDeintervalise(datalong=ndlong))
    if(timestep !='asdata') ndlong <- ctDiscretiseData(dlong=ndlong,timestep=timestep,
      TDpredNames=gm$TDpredNames,TIpredNames=gm$TIpredNames)
    
  if(gm$n.TDpred > 0) {
    gm$TDPREDMEANS=matrix(ndlong[,gm$TDpredNames],ncol=1)
    gm$TDPREDVAR = diag(0,gm$Tpoints*gm$n.TDpred)
    # gm$T0TDPREDCOV = matrix(0,gm$n.latent,gm$n.TDpred*gm$Tpoints)
    # gm$TRAITTDPREDCOV = matrix(0,gm$n.latent,gm$n.TDpred*gm$Tpoints)
  }
  if(gm$n.TIpred > 0) {
    gm$TIPREDMEANS=matrix(ndlong[,gm$TIpredNames],ncol=1)
    gm$TIPREDVAR = diag(0,gm$n.TIpred)
  }
  }
  if(timestep=='asdata') dtmat <- dat[,paste0('dT',1:(fit$ctmodelobj$Tpoints-1)),drop=FALSE] else dtmat <- NA
  new=suppressMessages(ctGenerate(ctmodelobj = gm,n.subjects = 1,dtmean=timestep,dtmat=dtmat,...))
  # new[,'id']=i
   out=rbind(out,new)
   # if(i==1 & n.subjects > 1) out=rbind(out,matrix(NA,nrow=nrow(out)*(n.subjects-1),ncol=ncol(out))) #preallocate
}
return(out)
}


