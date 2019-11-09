
#' Compares model implied density and values to observed, for a ctStanFit object.
#'
#' @param fit ctStanFit object. 
#' @param legend Logical, whether to plot a legend.
#' @param diffsize Integer > 0. Number of discrete time lags to use for data viz.
#' @param shading Logical -- show smoothed shading over generated data points? Otherwise, plot shaded polygon based on quantile estimate. 
#' Shading is better for non-linearities.
#' @param probs Vector of length 3 containing quantiles to plot -- should be rising numeric values between 0 and 1. 
#' @param wait Logical, if TRUE, waits for input before plotting next plot.
#' @param jitter Positive numeric between 0 and 1, if TRUE, jitters empirical data by specified proportion of std dev.
#' @param datarows integer vector specifying rows of data to plot. Otherwise 'all' uses all data.
#' @param ... extra arguments to pass to plot function.
#' @return If plot=FALSE, an array containing quantiles of generated data. If plot=TRUE, nothing, only plots.
#' @export
#' @details This function relies on the data generated during each iteration of fitting to approximate the
#' model implied distributions -- thus, when limited iterations are available, the approximation will be worse.
#'
#' @examples
#' \dontrun{
#' ctStanPostPredict(ctstantestfit,wait=FALSE, shading=FALSE, datarows=1:25,diffsize=2)
#' }
ctStanPostPredict <- function(fit,legend=TRUE,diffsize=1,jitter=.02, wait=TRUE,probs=c(.025,.5,.975),shading=TRUE, datarows='all',...){
 
  if(datarows[1]=='all') datarows <- 1:nrow(fit$data$Y)

  xmeasure=datarows
  if(is.null(fit$generated$Y)) Ygen <- ctStanGenerateData(fit)$generated$Y else Ygen <- fit$generated$Y
  Ygen<-aperm(Ygen,c(2,1,3))
  Ygen <- Ygen[,datarows,,drop=FALSE]
  time <- fit$standata$time[datarows]
  Ydat <- fit$data$Y[datarows,,drop=FALSE]
  ctDensityList(x=list(Ydat[datarows,,drop=FALSE],Ygen[,,,drop=FALSE]),plot=TRUE,
    main='All variables',lwd=2,legend = c('Observed','Model implied'),xlab='Value',...)
  
  y<-aaply(Ygen,c(2,3),quantile,na.rm=TRUE,probs=probs,.drop=FALSE)
  # y<-array(y,dim=dim(y)[-1])  
  dimnames(y) <- list(NULL,fit$ctstanmodel$manifestNames,paste0(probs*100,'%'))
  
  ycut=quantile(Ygen,c(.005,.995),na.rm=TRUE)
  ys=Ygen
  ys[ys<ycut[1] | ys>ycut[2]] <- NA 
  
  for(typei in c('obs','change')){
    for(i in 1:dim(Ygen)[3]){ #dim 3 is indicator, dim 2 is datarows, dim 1 iterations
      if(wait) readline("Press [return] for next plot.")
      # xsmeasure=rep(xmeasure,each=dim(Ygen)[1])
      # xstime=rep(xtime,each=dim(Ygen)[1])
      # xsmeasure=xsmeasure[ys>ycut[1] & ys<ycut[2]]
      # xstime=xstime[ys>ycut[1] & ys<ycut[2]]
      # 
      
      if(typei=='obs'){
        
        ctDensityList(x=list(Ydat[,i,drop=FALSE],Ygen[,datarows,i,drop=FALSE]),plot=TRUE,
          main=fit$ctstanmodel$manifestNames[i],lwd=2,legend = c('Observed','Model implied'),xlab='Value',...)
        
        if(wait) readline("Press [return] for next plot.")
        
        for(subtypei in c('Time','Observation')){
          if(subtypei=='Observation') x <- xmeasure
          if(subtypei=='Time') x <- time
          
          notmissing <- which(!is.na(c(y[datarows,i,1])))
          
          if(shading) {
            xs=rep(x,each=dim(Ygen)[1])
            ycut=quantile(Ygen[,,i],c(.005,.995),na.rm=TRUE)
            ysamps=Ygen[,,i]
            xs=xs[ysamps>ycut[1] & ysamps<ycut[2]]
            ysamps=ysamps[ysamps>ycut[1] & ysamps<ycut[2]]
            graphics::smoothScatter(xs,ysamps,nbin=256,colramp=grDevices::colorRampPalette(colors=c(rgb(1,1,1,0),rgb(1,.4,.4,.3))),nrpoints=0,
              transformation=function(x) x,ylim=range(c(y[,i,],quantile(ysamps,probs = c(.01,.99),na.rm=TRUE)),na.rm=TRUE),
              xlab=subtypei,ylab=dimnames(y)[[2]][i])
          }
          
          if(subtypei=='Observation') ctPlotArray(list(y=y[notmissing,i,,drop=FALSE],x=x[notmissing]),legend=FALSE,add=shading,polygon=!shading,
            plotcontrol=list(xlab=subtypei,main=dimnames(y)[[2]][i],...))
          
          # if(subtypei=='Time')
          
          ocol <- rgb(0,0,.7,.7)
          points(x[notmissing],
            Ydat[,i][notmissing] +  rnorm(length(Ydat[,i][notmissing]),0, jitter * sd(Ydat[,i][notmissing],na.rm=TRUE)),
            type=ifelse(subtypei=='Time','p','l'),lwd=2,lty=1,pch=17, col=ocol)
          if(legend) legend('topright',c('Model implied','Observed'),text.col=c('red',ocol))
          if(i < dim(Ygen)[3])  if(wait) readline("Press [return] for next plot.")
        }
      }
      
      
      if(typei=='change'){

        yp<-aperm(ys[,,i,drop=TRUE],c(2,1))#drop true set here if looking for problems!
        
        for(cdiffsize in diffsize){
          # diffindex <- c() 
          # for(diffi in 1:cdiffsize){
          #   diffindex <- c(diffindex,which(as.logical(fit$data$T0check[datarows][-1]))-(diffi-1))
          # }
          
          subdiff <- c( which(diff(fit$data$subject[datarows],cdiffsize) != 0), datarows[length(datarows):1][1:cdiffsize])
          
          dygen<-diff(yp,lag = cdiffsize) 
          # yp[-1,i,,,drop=FALSE] - yp[-fit$data$ndatapoints,i,,,drop=FALSE]
          dygendt <- dygen / diff(time,lag = cdiffsize)
          dygendt<-dygendt[-subdiff,,drop=FALSE]

          # dydt<-diff(Ydat[,i], lag = cdiffsize)/diff(time,lag = cdiffsize)
          dydt <- diff(Ydat[,i],lag=cdiffsize)
          dydt <- dydt[-subdiff]
          dydt <- dydt +  rnorm(length(dydt),0, jitter * sd(Ydat[,i],na.rm=TRUE)) #add jitter
          xtime <- time[-subdiff]
          # smoothScatter(matrix(yp[-fit$data$ndatapoints,i,,,drop=FALSE],ncol=1),
          #   matrix(dygendt[,,,,drop=FALSE],ncol=1),
          #   nbin=512,colramp=colorRampPalette(colors=c('white',rgb(1,0,0,1))),nrpoints=0,
          #   # transformation=function(x) x,
          #   ylim=range(c(quantile(c(dygendt),probs = c(.01,.99),na.rm=TRUE)),na.rm=TRUE),
          #   xlab='Observation',ylab=dimnames(y)[[2]][i])
          samps<-sample(1:length(dygendt),size=50000,replace=TRUE)
          plot(matrix(yp[-subdiff,,drop=FALSE][samps],ncol=1),
            matrix(dygendt[,,drop=FALSE][samps],ncol=1),
            ylab=paste0('dy/dt, diff=',cdiffsize),xlab='y', main=dimnames(y)[[2]][i],
            pch=16,cex=.2,col=rgb(1,0,0,.1),...)
          points( Ydat[-subdiff,i],
            dydt,
            col=rgb(0,0,1,.5),pch=17,...)
          
          if(wait) readline("Press [return] for next plot.")  
          
          plot(
            rep(xtime,(dim(dygendt)[2]))[samps],
            matrix(dygendt[,,drop=FALSE],ncol=1)[samps],
            ylab=paste0('dy/dt, diff=',cdiffsize),xlab='time', main=dimnames(y)[[2]][i],
            pch=16,cex=.1,col=rgb(1,0,0,.3),...)
          points(xtime,
            dydt,
            col=rgb(0,0,1,.5),pch=17,...)
          
        }
      }
    }
  }
  
}
