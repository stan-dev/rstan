#' summary.ctStanFit
#'
#' Summarise a ctStanFit object that was fit using \code{\link{ctStanFit}}. 
#' 
#' @param object fit object from \code{\link{ctStanFit}}, of class ctStanFit.
#' @param timeinterval positive numeric indicating time interval to use for discrete time parameter calculations
#' reported in summary. 
#' @param digits integer denoting number of digits to report.
#' @param parmatrices if TRUE, also return additional parameter matrices -- can be slow to compute
#' for large models with many samples.
#' @param priorcheck Whether or not to use \code{ctsem:::priorchecking} to compare posterior mean and sd to prior mean and sd.
#' @param ... Additional arguments to pass to \code{ctsem:::priorcheckreport}, such as \code{meanlim}, or \code{sdlim}.
#' @return List containing summary items.
#' @examples
#' summary(ctstantestfit)
#' @method summary ctStanFit
#' @export

summary.ctStanFit<-function(object,timeinterval=1,digits=3,parmatrices=FALSE,priorcheck=TRUE,...){
  
  if(class(object) != 'ctStanFit') stop('Not a ctStanFit object!')
  
  out=list()
  monvars <- c('mean','sd','2.5%','50%','97.5%')
  
  if(class(object$stanfit)=='stanfit'){ 
    s<-suppressWarnings(getMethod('summary','stanfit')(object$stanfit))
    if('98%' %in% colnames(s$summary)) colnames(s$summary)[colnames(s$summary)=='98%'] <- '97.5%'
    e <- extract(object) 
  }
  
  if(class(object$stanfit)!='stanfit')  e <- extract(object) 
  parnames <- c()
  parindices <- c()
  # for(m in names(object$setup$matsetup)){
  #   if(dim(object$setup$matsetup[[m]])[1] > 0){
  #     parnames <- c(parnames,rownames(object$setup$matsetup[[m]]))
  #     parindices <- c(parindices, object$setup$matsetup[[m]][,'param'])
  #   }
  # }

  parnames <- object$setup$popsetup$parname
  parindices <- object$setup$popsetup$param
  pars <- cbind(parnames,parindices)
  pars<-pars[!duplicated(pars[,1,drop=FALSE]),,drop=FALSE]
  parnames <- pars[as.numeric(pars[,2,drop=FALSE]) >0, 1]
  # parnames <- unique(parnames)
  parnamesiv <- parnames[object$data$indvaryingindex]
  
  #### generate covcor matrices of raw and transformed subject level params
  
  iter=dim(e$rawpopcorr)[1]
  if(!is.null(iter)){ #then there is some individual variation so continue
    nindvarying=dim(e$rawpopcorr)[2]
    
    if(nindvarying>1){
      
      getMean=function(myarray){
        out=matrix(NA,nrow=nindvarying,ncol=nindvarying)
        for(i in 1:nrow(out)){
          for(j in 1:ncol(out)){
            out[i,j]<-mean(myarray[i,j,])
          }}
        return(out)
      }
      
      getSd=function(myarray){
        out=matrix(NA,nrow=nindvarying,ncol=nindvarying)
        for(i in 1:nrow(out)){
          for(j in 1:ncol(out)){
            out[i,j]<-sd(myarray[i,j,])
          }}
        return(out)
      }
      
      if(1==99 & (is.null(object$data$intoverpop) || object$data$intoverpop==0)){    
        #transformed subject level params
        rawpopcorr_transformed= array(sapply(1:iter, function(x) cor(e$indparams[x,,])),dim=c(nindvarying,nindvarying,iter))
        rawpopcov_transformed= array(sapply(1:iter, function(x) cov(e$indparams[x,,])),dim=c(nindvarying,nindvarying,iter))
        
        rawpopcorr_transformedmean=getMean(rawpopcorr_transformed)
        rawpopcorr_transformedsd=getSd(rawpopcorr_transformed)
        
        rawpopcov_transformedmean=getMean(rawpopcov_transformed)
        rawpopcov_transformedsd=getSd(rawpopcov_transformed)
        
        rawpopcovcor_transformedmean=rawpopcov_transformedmean
        rawpopcovcor_transformedmean[lower.tri(diag(nindvarying))]=rawpopcorr_transformedmean[lower.tri(diag(nindvarying))]
        
        rawpopcovcor_transformedsd=rawpopcov_transformedsd
        rawpopcovcor_transformedsd[lower.tri(diag(nindvarying))]=rawpopcorr_transformedsd[lower.tri(diag(nindvarying))]
        
        dimnames(rawpopcovcor_transformedsd)<-list(parnamesiv,parnamesiv)
        dimnames(rawpopcovcor_transformedmean)<-list(parnamesiv,parnamesiv)
        
        out=list(paste0('The following matrix is the posterior mean of the correlation and covariance matrix of subject level parameters,', 
          ' with correlations on the lower triangle'),
          popcovcor_mean=round(rawpopcovcor_transformedmean,digits),
          paste('The following matrix is the posterior std dev. of the correlation and covariance matrix of subject level parameters,', 
            'with correlations on the lower triangle'),
          popcovcor_sd=round(rawpopcovcor_transformedsd,digits))
      }
      #raw pop distribution params
      dimrawpopcorr <- dim(e$rawpopcorr)
      if(class(object$stanfit)!='stanfit') rawpopcorr= array(e$rawpopcorr,dim=c(dimrawpopcorr[1],1,dimrawpopcorr[2] * dimrawpopcorr[3]))
      if(class(object$stanfit)=='stanfit') rawpopcorr= rstan::extract(object$stanfit,pars='rawpopcorr',permuted=FALSE)
      
      rawpopcorrout <- suppressWarnings(monitor(rawpopcorr, digits_summary=digits,warmup=0,print = FALSE)[lower.tri(diag(nindvarying)),c(monvars,'n_eff','Rhat'),drop=FALSE])
      
      # rawpopcorrout <- ctCollapse(rawpopcorr,1,mean)
      # rawpopcorrout <- cbind(rawpopcorrout,ctCollapse(rawpopcorr,1,sd)[lower.tri(diag(nindvarying)),drop=FALSE])
      # rawpopcorrout <- cbind(rawpopcorrout,ctCollapse(rawpopcorr,1,quantile,probs=c(.025))[lower.tri(diag(nindvarying)),drop=FALSE])
      # rawpopcorrout <- cbind(rawpopcorrout,ctCollapse(rawpopcorr,1,quantile,probs=c(.5))[lower.tri(diag(nindvarying)),drop=FALSE])
      # rawpopcorrout <- cbind(rawpopcorrout,ctCollapse(rawpopcorr,1,quantile,probs=c(.975))[lower.tri(diag(nindvarying)),drop=FALSE])
      # colnames(rawpopcorrout) <- monvars
      rownames(rawpopcorrout) <- matrix(paste0('',parnamesiv,'__',rep(parnamesiv,each=length(parnamesiv))),
        length(parnamesiv),length(parnamesiv))[lower.tri(diag(nindvarying)),drop=FALSE]
      # rawpopcorrout <- round(rawpopcorrout,digits=digits)
      
      rawpopcorrout <- cbind(rawpopcorrout,rawpopcorrout[,'mean'] / rawpopcorrout[,'sd'])
      colnames(rawpopcorrout)[ncol(rawpopcorrout)] <- 'z'
      
      # rawpopcorrout <- rawpopcorrout[order(abs(rawpopcorrout[,'z'])),,drop=FALSE]
      
      rawpopcorrmean= ctCollapse(e$rawpopcorr,1,mean)
      rawpopcorrsd= ctCollapse(e$rawpopcorr,1,sd)
      rawpopcov_mean = ctCollapse(e$rawpopcov,1,mean)
      rawpopcov_sd=ctCollapse(e$rawpopcov,1,sd)
      
      dimnames(rawpopcorrmean)<-list(parnamesiv,parnamesiv)
      dimnames(rawpopcorrsd)<-list(parnamesiv,parnamesiv)
      dimnames(rawpopcov_mean)<-list(parnamesiv,parnamesiv)
      dimnames(rawpopcov_sd)<-list(parnamesiv,parnamesiv)
      
      out=list(note='Posterior means of the raw parameter population distribution correlation matrix:',
        rawpopcorr_mean=round(rawpopcorrmean,digits),
        note='Posterior std dev. of the raw parameter population distribution correlation matrix:',
        rawpopcorr_sd=round(rawpopcorrsd,digits),
        note='Posterior means of the raw parameter population distribution covariance matrix:',
        rawpopcov_mean = round(rawpopcov_mean,digits),
        note='Posterior std dev. of the raw parameter population distribution covariance matrix:',
        rawpopcov_sd=round(rawpopcov_sd,digits)
      )
      
      out$rawpopcorr = round(rawpopcorrout,digits)
    }
  }
  
  if(priorcheck & object$standata$nopriors==0) out = c(out,priorcheckreport(object,...))
  
  if(object$ctstanmodel$n.TIpred > 0) {

    if(class(object$stanfit)=='stanfit'){
      rawtieffect <- rstan::extract(object$stanfit,permuted=FALSE,pars='TIPREDEFFECT')
      tidiags <- suppressWarnings(monitor(rawtieffect,warmup=0,digits_summary = digits,print = FALSE))
    }
    tieffect <- array(e$linearTIPREDEFFECT,dim=c(dim(e$linearTIPREDEFFECT)[1], 1, length(parnames) * dim(e$linearTIPREDEFFECT)[3]))
    tieffectnames <- paste0('tip_',rep(object$ctstanmodel$TIpredNames,each=length(parnames)),'_',parnames)
    dimnames(tieffect)<-list(c(),c(),tieffectnames)
    tipreds = suppressWarnings(monitor(tieffect,warmup = 0,print = FALSE)[,monvars])
    if(class(object$stanfit)=='stanfit') tipreds <- cbind(tipreds,tidiags[,c('n_eff','Rhat')])
    tipreds <- tipreds[c(object$data$TIPREDEFFECTsetup)>0,,drop=FALSE]
    z = tipreds[,'mean'] / tipreds[,'sd'] 
    out$tipreds= round(cbind(tipreds,z),digits) #[order(abs(z)),]
  }
  
  if(parmatrices){
    
    # #check if stanfit object can be used
    # sf <- object$stanfit
    # npars <- try(get_num_upars(sf),silent=TRUE) #$stanmodel)
    # 
    # if(class(npars)=='try-error'){ #in case R has been restarted or similar
    #   standataout <- object$standata
    #   smf <- stan_reinitsf(object$stanmodel,standataout)
    # }
  object$standata$savescores <- 0L
  object$standata$gendata <- 0L
  object$standata$dokalman <- 0L
  sf <- stan_reinitsf(object$stanmodel,data=object$standata)

    parmatlists <- try(apply(e$rawpopmeans[
      sample(x = 1:dim(e$rawpopmeans)[1],
        size = dim(e$rawpopmeans)[1],  #min(parmatsamples,
        replace = FALSE),],
      1,ctStanParMatrices,fit=object,timeinterval=timeinterval,sf=sf))
      
    if(class(parmatlists)!='try-error'){
      parmatarray <- array(unlist(parmatlists),dim=c(length(unlist(parmatlists[[1]])),length(parmatlists)))
      parmats <- matrix(NA,nrow=length(unlist(parmatlists[[1]])),ncol=7)
      rownames(parmats) <- paste0('r',1:nrow(parmats))
      counter=0
      for(mati in 1:length(parmatlists[[1]])){
        if(all(dim(parmatlists[[1]][[mati]]) > 0)){
        for(coli in 1:ncol(parmatlists[[1]][[mati]])){
          for(rowi in 1:nrow(parmatlists[[1]][[mati]])){
            counter=counter+1
            new <- matrix(c(
              rowi,
              coli,
              mean(parmatarray[counter,],na.rm=TRUE),
              sd(parmatarray[counter,],na.rm=TRUE),
              quantile(parmatarray[counter,],probs=c(.025,.5,.975),na.rm=TRUE)),
              nrow=1)
            try(rownames(parmats)[counter] <- names(parmatlists[[1]])[mati])
            try(parmats[counter,]<-new)
          }}}}
      colnames(parmats) <- c('Row','Col', 'Mean','Sd','2.5%','50%','97.5%')
      
      #remove certain parmatrices lines
      removeindices <- which(rownames(parmats) == 'MANIFESTVAR' & parmats[,'Row'] != parmats[,'Col'])
      
      removeindices <- c(removeindices,which((rownames(parmats) %in% c('MANIFESTVAR','T0VAR','DIFFUSION','dtDIFFUSION','asymDIFFUSION',
        'T0VARcor','DIFFUSIONcor','dtDIFFUSIONcor','asymDIFFUSIONcor') &  parmats[,'Row'] < parmats[,'Col'])))
      
      removeindices <- c(removeindices,which((rownames(parmats) %in% c('T0VARcor','DIFFUSIONcor','dtDIFFUSIONcor','asymDIFFUSIONcor') & 
          parmats[,'Row'] == parmats[,'Col'])))
      
      parmats <- parmats[-removeindices,]
      
      
      out$parmatrices=round(parmats,digits=digits)
      
      out$parmatNote=paste0('Population mean parameter matrices calculated with time interval of ', timeinterval,' for discrete time (dt) matrices. ',
        'Covariance related matrices shown as covariance matrices, correlations have (cor) suffix. Asymptotic (asym) matrices based on infinitely large time interval.')
    }
    if(class(parmatlists)=='try-error') out$parmatNote = 'Could not calculate parameter matrices'
  }
  
  
  if(class(object$stanfit)=='stanfit'){
    popsd=s$summary[c(grep('^popsd',rownames(s$summary),fixed=FALSE)),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE] [ object$data$indvaryingindex,,drop=FALSE]
    rownames(popsd)=parnames[ object$data$indvaryingindex]
    # popmeans=s$summary[c(grep('hmean_',rownames(s$summary))),
    #   c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE]
    popmeans=s$summary[c(grep('popmeans[', rownames(s$summary),fixed=TRUE)),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE]
    popmeans=popmeans[(nrow(popmeans)/2+1):nrow(popmeans),,drop=FALSE]
    rownames(popmeans) <- parnames
    
    
    
    logprob=s$summary[c(grep('lp',rownames(s$summary))),
      c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE]
  }
  
  if(class(object$stanfit)!='stanfit'){ #if optimized / importance sampled

    if(!is.null(iter)){ popsd <- suppressWarnings(monitor(array(e$popsd,dim=c(dim(e$popsd)[1],1,dim(e$popsd)[2])),warmup=0,print=FALSE))
    popsd=popsd[ object$data$indvaryingindex, monvars,drop=FALSE]
    rownames(popsd)=parnamesiv
    }

    popmeans=suppressWarnings(monitor(array(e$popmeans,dim=c(dim(e$popmeans)[1],1,dim(e$popmeans)[2])),warmup=0,print=FALSE))
    rownames(popmeans) = parnames #names(e)[grep('hmean_',names(e))]
    popmeans = popmeans[,monvars,drop=FALSE]
    
    logprob = object$stanfit$optimfit$value
    aic = 2*ncol(object$stanfit$rawposterior) - 2*logprob
  }
  
  if(!is.null(iter)) out$popsd=round(popsd,digits=digits)
  
  out$popmeans=round(popmeans,digits=digits)
  
  out$popNote=paste0('popmeans are reported as specified in ctModel -- covariance related matrices are in sd / matrix square root form.')
  
  out$logprob=logprob
  
  if(class(object$stanfit)!='stanfit') out$aic = aic
  
  if(!parmatrices) out$parmatNote <- 'For additional summary matrices, use argument: parmatrices = TRUE'
  
  
  
  # out$posteriorpredictive=round(s$summary[c(grep('stateppll',rownames(s$summary))),
  #     c('mean','sd','2.5%','50%','97.5%','n_eff','Rhat'),drop=FALSE],3)
  # }
  
  
  # if(class(object$stanfit)!='stanfit'){ #optimization summary
  #   out=list()
  #   out$popmeans=object$stanfit$transformedpars[grep('hmean_',rownames(object$stanfit$transformedpars)),]
  #   out$popsd=object$stanfit$transformedpars[grep('hsd_',rownames(object$stanfit$transformedpars)),]
  #   out$logprob=round(-object$stanfit$optimfit$value)
  # }
  
  return(out)
}
