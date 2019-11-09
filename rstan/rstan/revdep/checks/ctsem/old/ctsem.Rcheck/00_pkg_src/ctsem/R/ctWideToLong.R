#' ctWideToLong
#' Convert ctsem wide to long format
#' @param datawide ctsem wide format data
#' @param Tpoints number of measurement occasions in data
#' @param n.manifest number of manifest variables
#' @param n.TDpred number of time dependent predictors
#' @param n.TIpred number of time independent predictors
#' @param manifestNames Character vector of manifest variable names.
#' @param TDpredNames Character vector of time dependent predictor names.
#' @param TIpredNames Character vector of time independent predictor names.
#' @details
#' Names must account for *all* the columns in the data - i.e. do not leave certain variables out
#' just because you do not need them.
#' @examples 
#'  #First load the example ctsem wide format data with absolute times
#'  data('datastructure')
#'  datastructure #contains two time intervals (dTx), therefore 3 time points.
#'  #Then convert to long format
#'  longexample <- ctWideToLong(datawide = datastructure, Tpoints=3, 
#'  n.manifest=3, manifestNames = c("Y1", "Y2", "Y3"),
#'  n.TDpred=1, TDpredNames = "TD1", 
#'  n.TIpred=2, TIpredNames = c("TI1", "TI2"))
#'
#'  #Then convert the time intervals to absolute time
#'  long <- ctDeintervalise(datalong = longexample, id='id', dT='dT')
#'  long
#'
#' 
#' @export

ctWideToLong<-function(datawide,Tpoints,n.manifest,n.TDpred=0, n.TIpred=0, 
  manifestNames='auto',TDpredNames='auto',TIpredNames='auto'){ 
  
  #names
  if(all(manifestNames=='auto')) manifestNames=paste0('Y',1:n.manifest)
  if(length(manifestNames) != n.manifest) stop("Length of manifestNames does not equal n.manifest!") 

  if(n.TDpred > 0){
    if(all(TDpredNames=='auto')) TDpredNames=paste0('TD',1:n.TDpred)
    if(length(TDpredNames) != n.TDpred) stop("Length of TDpredNames does not equal n.TDpred!") 
  }
  
  if(n.TIpred > 0){
    if(all(TIpredNames=='auto')) TIpredNames=paste0('TI',1:n.TIpred)
    if(length(TIpredNames) != n.TIpred) stop("Length of TIpredNames does not equal n.TIpred!") 
  }
  
  datawide<-as.matrix(datawide,nrow=nrow(datawide),ncol=ncol(datawide)) #set to matrix for manipulations
  n.subjects<-nrow(datawide) #calculate number of subjects in datawide
  
  
  manifests<-matrix(t(datawide[,
    paste0(manifestNames,'_T',rep(0:(Tpoints-1),each=n.manifest)),
      drop=FALSE]),byrow=T,nrow=n.subjects*Tpoints)
  colnames(manifests)<-manifestNames
  
  times<-matrix(t(cbind(0,datawide[,(Tpoints*n.manifest+(Tpoints)*n.TDpred+1) : 
      (Tpoints*n.manifest+(Tpoints)*n.TDpred+ (Tpoints-1)),drop=FALSE])),byrow=T,nrow=n.subjects*Tpoints)
  colnames(times)<-'dT'
  
  id<-matrix(rep(1:n.subjects,each=Tpoints),ncol=1)
  colnames(id)<-'id'
  
  datalong<-cbind(id,times,manifests)
  
  if(n.TDpred>0) {
    
    tdpreds<-matrix(t(datawide[,
      paste0(TDpredNames,'_T',rep(0:(Tpoints-1),each=n.TDpred)),
      drop=FALSE]),byrow=T,nrow=n.subjects*Tpoints)
    colnames(tdpreds)<-TDpredNames
    datalong<-cbind(datalong,tdpreds)
    
  # tdpreds<-matrix(NA,nrow=n.subjects*(Tpoints),ncol=n.TDpred)
  # for(tdpredi in 1:n.TDpred){
  # tdpreds[,tdpredi]<-t(datawide[,(Tpoints*n.manifest+(Tpoints)*(tdpredi-1)+1) : (Tpoints*n.manifest+(Tpoints)*(tdpredi)),drop=FALSE])
  # }
  # tdpredsfull<-matrix(NA,nrow=n.subjects*Tpoints,ncol=n.TDpred)
  # tdpredsfull[,]<-tdpreds
  #   colnames(tdpredsfull)<-TDpredNames
  #   datalong<-cbind(datalong,tdpredsfull)
    
  }
  
  if(n.TIpred > 0){
   
  tipreds<- matrix(t(datawide[,(Tpoints*n.manifest+(Tpoints)*n.TDpred+(Tpoints)) : 
      (Tpoints*n.manifest+(Tpoints)*n.TDpred+ Tpoints -1 + n.TIpred),drop=FALSE]),byrow=T,nrow=n.subjects)
  tipredsfull<-matrix(NA,nrow=n.subjects*Tpoints,ncol=n.TIpred)
  tipredsfull<-tipreds[rep(1:n.subjects,each=Tpoints),,drop=F]
  colnames(tipredsfull)<-TIpredNames
  datalong<-cbind(datalong,tipredsfull)
  }
    
  return(datalong)
}
