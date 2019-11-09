#' Get time independent predictor effect estimates
#' 
#' Computes and plots combined effects and quantiles for effects of time independent predictors
#' on subject level parameters of a ctStanFit object.
#'
#' @param fit fit object from \code{\link{ctStanFit}}
#' @param returndifference logical. If FALSE, absolute parameter values are returned. 
#' If TRUE, only the effect of the covariate (i.e. without the average value of the parameter)
#' are returned. The former can be easier to interpret, but the latter are more likely to fit multiple plots together. 
#' Not used if \code{parmatrices=TRUE}.
#' @param probs numeric vector of quantile probabilities from 0 to 1. Specify 3
#' values if plotting, the 2nd will be drawn as a line with uncertainty polygon
#' based on 1st and 3rd.
#' @param includeMeanUncertainty if TRUE, output includes sampling variation in the mean parameters. If FALSE,
#' mean parameters are fixed at their median, only uncertainty in time independent predictor effects is included.  
#' @param whichTIpreds integer vector specifying which of the tipreds in the fit object you want to
#' use to calculate effects. Unless quadratic / higher order versions of predictors have been 
#' included, selecting more than one probably doesn't make sense. If for instance a squared
#' predictor has been included, then you can specify both the linear and squared version. 
#' The x axis of the plot (if generated) will be based off the first indexed predictor. To 
#' check what predictors are in the model, run \code{fit$ctstanmodel$TIpredNames}.
#' @param parmatrices Logical. If TRUE (default), the \code{\link{ctStanParMatrices}} function
#' is used to return an expanded range of possible matrices of interest.
#' @param whichpars if parmatrices==TRUE, character vector specifying which matrices, and potentially which 
#' indices of the matrices, to plot. c('dtDRIFT[2,1]', 'DRIFT') would output for row 2 and column 1 of 
#' the discrete time drift matrix, as well as all indices of the continuous time drift matrix. 
#' If parmatrices==FALSE, integer vector specifying which of the subject
#' level parameters to compute effects on. The integers corresponding to certain parameters can be found in the 
#' \code{param} column of the \code{fit$setup$popsetup} object. In either case 'all' uses all available parameters.
#' @param nsamples Positive integer specifying the maximum number of saved iterations to use. 
#' Character string 'all' can also be used.
#' @param nsubjects Positive integer specifying the number of subjects to compute values for.
#' Character string 'all' can also be used. Time taken is a function of nsubjects*niterations.
#' @param plot Logical. If TRUE, nothing is returned but instead \code{\link{ctPlotArray}}
#' is used to plot the output instead.
#' @param timeinterval positive numeric indicating time interval to use for discrete time parameter matrices,
#' if \code{parmatrices=TRUE}.
#' @param filter either NA, or a length 2 vector, where the first element contains the time independent predictor index
#' to filter by, and the second contains the comparison operator in string form (e.g. "< 3",
#' to only calculate effects for subjects where the tipreds of the denoted index are less than 3).
#' @param ... arguments to pass to \code{\link{ctPlotArray}} for plotting.
#' @return Either a three dimensional array of predictor effects, or nothing with a plot
#' generated.
#' @export
#'
#' @examples
#' \dontrun{
#' #samples reduced here for speed
#' ctStanTIpredeffects(ctstantestfit,plot=TRUE,whichpars='CINT',nsamples=10,nsubjects=10)
#' }
ctStanTIpredeffects<-function(fit,returndifference=FALSE, probs=c(.025,.5,.975),
  includeMeanUncertainty=FALSE,
  whichTIpreds=1,parmatrices=TRUE, whichpars='all', nsamples=100, timeinterval=1,
  nsubjects=50,filter=NA,
  plot=FALSE,...){

  ctspec <- fit$ctstanmodel$pars
  
  e<-extract(fit)
  rawpopmeans <- e$rawpopmeans
  
  niter<-dim(e$rawpopmeans)[1]
  if(nsamples=='all' || nsamples > niter) nsamples <- niter
  rawpopmeans <- rawpopmeans[sample(x = 1:niter, nsamples,replace = FALSE),]
  
  message(sprintf('Getting %s samples by %s subjects for %s total samples', nsamples,nsubjects, nsamples*nsubjects))
  
  if(!includeMeanUncertainty) rawpopmeans <- matrix(apply(rawpopmeans,2,median),byrow=TRUE,nrow=nrow(rawpopmeans),ncol=ncol(rawpopmeans))
  
  tipreds<-ctCollapse(e$tipreds,1,mean) #maybe collapsing over sampled tipred values is not ideal?
  
  if(any(!is.na(filter))) tipreds <- eval(parse(text=paste0('tipreds[tipreds[,',filter[1],']',filter[2],',,drop=FALSE]')))
  tipreds <- tipreds[,whichTIpreds,drop=FALSE]
  if(nsubjects=='all') nsubjects = nrow(tipreds)
  if(nsubjects > nrow(tipreds) && length(whichTIpreds)>1) nsubjects <- nrow(tipreds) #if we need to use real tipred data, no point using more subjects
  
  #use tipred data or generated points? latter is smoother but can't use in case of interactions.
  if(length(whichTIpreds)>1) tipreds <- tipreds[sample(x = 1:nsubjects, nsubjects,replace = FALSE),,drop=FALSE]
  if(length(whichTIpreds)==1) tipreds<-cbind(seq(from=min(tipreds), to = max(tipreds,na.rm=TRUE), length.out=nsubjects))
  
  
  tieffect<-e$TIPREDEFFECT[,,whichTIpreds,drop=FALSE]

  tiorder<-order(tipreds[,1])
  tipreds<-tipreds[tiorder,,drop=FALSE] #order tipreds according to first one
  
  message('Calculating time independent predictor effects...')
  
  raweffect <- aaply(1:nrow(rawpopmeans),1,function(iterx) { #for every iter
    aaply(tipreds,1,function(tix){ #and every distinct tipred vector
      rawpopmeans[iterx,,drop=FALSE] + t(matrix(tieffect[iterx,,,drop=FALSE],nrow=dim(tieffect)[2]) %*% tix)
    },.drop=FALSE)
  },.drop=FALSE)
  
  
  if(!parmatrices) {
    if(all(whichpars=='all')) whichpars=which(apply(tieffect,2,function(x) any(x!=0)))
    raweffect <- raweffect[,,,whichpars,drop=FALSE]
    tieffect<-tieffect[,whichpars,,drop=FALSE] #updating...
    npars<-length(whichpars)
    effect<-aaply(1:npars, 1,function(pari){ #for each param
      param=raweffect[,,,pari]
      out=tform(param, 
        fit$setup$popsetup$transform[whichpars[pari]], 
        fit$setup$popvalues$multiplier[whichpars[pari]],
        fit$setup$popvalues$meanscale[whichpars[pari]],
        fit$setup$popvalues$offset[whichpars[pari]], 
        fit$setup$popvalues$inneroffset[whichpars[pari]], fit$setup$extratforms) 
      return(out)
    })
    if(returndifference){ #if only returning differences from zero
      noeffect<-aaply(1:npars, 1,function(pari){ #for each param
        param <- rawpopmeans[,pari]
        out=tform(param, 
        fit$setup$popsetup$transform[whichpars[pari]], 
        fit$setup$popvalues$multiplier[whichpars[pari]],
        fit$setup$popvalues$meanscale[whichpars[pari]],
        fit$setup$popvalues$offset[whichpars[pari]], 
        fit$setup$popvalues$inneroffset[whichpars[pari]], fit$setup$extratforms) 
        return(out)
      })
      effect<-effect-array(noeffect,dim=dim(effect))
    }
  }
  
  
  if(parmatrices)  {
    rawpopmeans <- rawpopmeans[rep(1:nrow(rawpopmeans),each=nsubjects),] #match rows of rawpopmeans and raweffect
    raweffect <- matrix(raweffect,ncol=dim(raweffect)[4])
    
    # #check if stanfit object can be used
    # sf <- fit$stanfit
    # npars <- try(get_num_upars(sf),silent=TRUE) #$stanmodel)
    # 
    # if(class(npars)=='try-error'){ #in case R has been restarted or similar
    #   sf <- stan_reinitsf(fit$stanmodel,fit$standata)
    # }

#adjust for speed
  fit$standata$savescores <- 0L
  fit$standata$gendata <- 0L
  fit$standata$dokalman <- 0L
  sf <- stan_reinitsf(fit$stanmodel,data=fit$standata)
  whichmatrices <- sapply(whichpars,function(x) {
    x=gsub(pattern = '[','',x,fixed=TRUE)
    x=gsub(pattern = ']','',x,fixed=TRUE)
    x=gsub(pattern = '[0-9]','',x)
    x=gsub(pattern = ',','',x,fixed=TRUE)
  })
    parmatlists<-lapply(1:nrow(rawpopmeans), function(x) { #for each param vector
      out = ctStanParMatrices(fit,  raweffect[x,], timeinterval=timeinterval,sf = sf)#rawpopmeans[x,] +
      return(out)
    })
    # browser()
    
    # parmatlists <- apply(e$rawpopmeans,1,ctStanParMatrices,model=object,timeinterval=timeinterval)
    parmatarray <- array(unlist(parmatlists),dim=c(length(unlist(parmatlists[[1]])),length(parmatlists)))
    parmats <- matrix(0,nrow=0,ncol=2)
    counter=0
    for(mati in 1:length(parmatlists[[1]])){
      if(all(dim(parmatlists[[1]][[mati]]) > 0)){
      for(coli in 1:ncol(parmatlists[[1]][[mati]])){
        for(rowi in 1:nrow(parmatlists[[1]][[mati]])){
          counter=counter+1
          new <- matrix(c(
            rowi,
            coli),
            nrow=1)
          rownames(new) = paste0(names(parmatlists[[1]])[[mati]])
          parmats<-rbind(parmats, new)
        }}}}
    colnames(parmats) <- c('Row','Col') 
    
    rownames(parmatarray) <- paste0(rownames(parmats),'[',parmats[,'Row'],',',parmats[,'Col'],']')
    
    #remove certain parmatrices lines
    removeindices <- which(rownames(parmats) == 'MANIFESTVAR' & parmats[,'Row'] != parmats[,'Col'])
    
    removeindices <- c(removeindices,which((rownames(parmats) %in% c('MANIFESTVAR','T0VAR','DIFFUSION','dtDIFFUSION','asymDIFFUSION',
      'T0VARcor','DIFFUSIONcor','dtDIFFUSIONcor','asymDIFFUSIONcor') &  parmats[,'Row'] < parmats[,'Col'])))
    
    removeindices <- c(removeindices,which((rownames(parmats) %in% c('T0VARcor','DIFFUSIONcor','dtDIFFUSIONcor','asymDIFFUSIONcor') & 
        parmats[,'Row'] == parmats[,'Col'])))
    
    parmatarray <- parmatarray[-removeindices,]
    
    effect <- array(parmatarray,dim=c(nrow(parmatarray),nsamples,nsubjects))
    rownames(effect) <- rownames(parmatarray)
    if(any(whichpars !='all')) {
      selection <- unlist(lapply(whichpars,function(x) grep(paste0('^\\Q',x,'\\E'),dimnames(effect)[[1]])))
      effect <- effect[selection,,,drop=FALSE]
    }
    
  }    
  
  
  out<-aaply(probs,1,function(x) ctCollapse(effect,2,quantile,probs=x,na.rm=TRUE),.drop=FALSE)
  
  if(!parmatrices) dimnames(out)=list(Quantile=paste0('Quantile',probs),
    popmean=fit$setup$popsetup$parname[whichpars],
    subject=tiorder #subjects reordered because tipreds were at top
  )
  if(parmatrices) dimnames(out)=list(Quantile=paste0('Quantile',probs),
    param=rownames(effect),
    subject=tiorder #subjects reordered because tipreds were at top
  )
  
  colnames(tipreds) <- colnames(fit$data$tipredsdata)[whichTIpreds]
  out <- list(y=aperm(out, c(3,2,1)), x=tipreds[,1,drop=FALSE])
  
  if(!plot) return(out) else {
    dots <- list(...)
    dots$input=out
    # dots$x=tipreds[,1]
    if(is.null(dots$plotcontrol)) dots$plotcontrol=list(
      ylab=ifelse(!returndifference,'Par. Value','Effect'),
      xlab=colnames(tipreds)[1],
      xaxs='i')
    
    do.call(ctPlotArray,dots)
  }
}

