#' ctCompareExpected
#' Compares model implied to observed means and covariances for panel data fit with ctsem. 
#' @param fitobj Fitted model object from OpenMx or ctsem.
#' @param cov Logical. If TRUE, show covariance plots, if FALSE show correlations.
#' @param outputmatrices if TRUE, output expected, observed, and residual correlation matrices as well as plots.
#' @param pause if TRUE (default) output plots interactively, one at a time.  If FALSE, output without stopping.
#' @param varlist if "all" include all variables in dataset.  Otherwise, specify numeric vector of variables to include.
#' @param ylim vector of min and max Y axis limits for plot.
#' @param ... additional arguments passed to plot. 
#' @export
ctCompareExpected<-function(fitobj,cov=TRUE,
  outputmatrices=FALSE,pause=TRUE,varlist="all",ylim=c(-1,1),...){ 

  if(class(fitobj)=="ctsemFit"){ 
    n.manifest<-fitobj$ctmodelobj$n.manifest
    Tpoints<-fitobj$ctmodelobj$Tpoints
    mxobj<-fitobj$mxobj    
  }
  
  if(class(fitobj)=="MxModel" &
      any( c(!exists('n.manifest'), !exists('Tpoints')))) message(
        'OpenMx object specified - additional arguments required: n.manifest, Tpoints') 
  

  if(varlist=="all") varlist <- 1:n.manifest
  
  datawide<-mxobj@data@observed
  if(dim(datawide)[1] == dim(datawide)[2]) original<-stats::cov2cor(datawide) else {
    datawide<-datawide[,paste0(fitobj$ctmodelobj$manifestNames, '_T',rep(0:(Tpoints-1),each=n.manifest))]
    if(cov!=TRUE) original<-round(stats::cor(datawide,use="pairwise.complete.obs"),digits=3) #subobtimal
    if(cov==TRUE) original<-round(stats::cov(datawide,use="pairwise.complete.obs"),digits=3) #subobtimal
  } 
  
if(cov==TRUE) modelexp<-mxobj$fitfunction$info$expCov[1:(n.manifest*Tpoints),1:(n.manifest*Tpoints)]
  if(cov!=TRUE) modelexp<-stats::cov2cor(mxobj$fitfunction$info$expCov[1:(n.manifest*Tpoints),1:(n.manifest*Tpoints)])
  
  expmean <- mxobj$fitfunction$info$expMean
  origmean <- apply(datawide,2,mean,na.rm=T)
  

  residuals <-  original-modelexp

  out<-list(original,modelexp,residuals)
  names(out)<-c("original","expected","residual")
  


  for(i in varlist){
    for(j in varlist){
      
      if(i==j) { #plot mean
      
        graphics::plot(0:(Tpoints-1),expmean[seq(i,Tpoints*n.manifest,n.manifest)],
          type='b',col='red',
          ylim=c(min(c(expmean[seq(i,Tpoints*n.manifest,n.manifest)], origmean[seq(i,Tpoints*n.manifest,n.manifest)]),na.rm=T),
            max(c(expmean[seq(i,Tpoints*n.manifest,n.manifest)], origmean[seq(i,Tpoints*n.manifest,n.manifest)]),na.rm=T)),
          ylab='Value',xlab='Measurement occasion', main=paste0('Means variable ',i))
        graphics::points(0:(Tpoints-1),origmean[seq(i,Tpoints*n.manifest,n.manifest)],col='blue',type='b')
        
        if(pause==T){
          message("Press [enter] to display next graph")
          readline()
        }
      }
    
      suppressWarnings(graphics::plot(rep(seq(n.manifest,Tpoints*n.manifest,n.manifest)/n.manifest,each=Tpoints)+ #x section
          seq(-.3,.3,length.out=Tpoints),#placement within x section
        (out$original[cbind(rep(seq(j,Tpoints*n.manifest,n.manifest),times=Tpoints),
          rep(seq(i,Tpoints*n.manifest,n.manifest),each=Tpoints))]),
        main=paste0(ifelse(cov==TRUE,'cov','cor')," variable ",i,"*",j),
        ylab="Value",xlab="Measurement occasion",col="blue",pch=16,
        ylim=c( ifelse(cov!=TRUE,-1,min(c(out$original,out$expected),na.rm=T)),
          ifelse(cov!=TRUE,1,max(c(out$original,out$expected),na.rm=T))),...))
      graphics::abline(v=seq(1.5,Tpoints-.5,1),col="grey",lty=2)
      
      
      graphics::points(rep(seq(n.manifest,Tpoints*n.manifest,n.manifest)/n.manifest,each=Tpoints)+
          seq(-.3,.3,length.out=Tpoints),
        out$expected[cbind(rep(seq(j,Tpoints*n.manifest,n.manifest),times=Tpoints),
          rep(seq(i,Tpoints*n.manifest,n.manifest),each=Tpoints))],col="red")
      
      graphics::legend("bottomleft",legend=list("original","expected"),text.col=c("blue","red"),bty="n")
      
      if(i*j<n.manifest^2){
        if(pause==T){
          message("Press [enter] to display next graph")
          readline()
        }
      }
    }
  }
  if(outputmatrices==TRUE) return(out)
}