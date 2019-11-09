#' Discretise long format continuous time (ctsem) data to specific timestep.
#'
#' Extends and rounds timing information so equal intervals, according to specified
#' timestep, are achieved. NA's are inserted in other columns as necessary,
#' any columns specified by TDpredNames or TIpredNames have zeroes rather than NA's
#' inserted (because some estimation routines do not tolerate NA's in covariates).
#'
#' @param dlong Long format data
#' @param timestep Positive real value to discretise
#' @param timecol Name of column containing absolute (not intervals) time information.
#' @param idcol Name of column containing subject id variable.
#' @param TDpredNames Vector of column names of any time dependent predictors
#' @param TIpredNames Vector of column names of any time independent predictors
#'
#' @return long format ctsem data.
#' @export
#'
#' @examples
#' long <- ctWideToLong(datawide=ctExample2,Tpoints=8,n.manifest=2,n.TDpred=1,
#'  manifestNames=c('LeisureTime','Happiness'),
#'  TDpredNames=c('MoneyInt'))
#' 
#' long <- ctDeintervalise(long)
#' 
#' long <- ctDiscretiseData(dlong=long, timestep = 1.1,TDpredNames=c('MoneyInt'))

ctDiscretiseData <- function(dlong,timestep,timecol='time',idcol='id',TDpredNames=NULL,
  TIpredNames=NULL){
  
  if(any(is.na(dlong[,timecol]))) stop('Cannot discretise with missing time data!')
  if(any(is.na(dlong[,idcol]))) stop('Cannot discretise with missing id data!')
  
  out<-matrix(NA,nrow=0,ncol=ncol(dlong))
  for(idi in unique(dlong[,idcol])){
    odat<-dlong[dlong[,idcol]==idi,]
    odat[,timecol]<-plyr::round_any(odat[,timecol],timestep)
    trange<-range(odat[,timecol])
    time<-seq(trange[1],trange[2],timestep)
    ndat<-matrix(NA,nrow=length(time),ncol=ncol(dlong),dimnames=dimnames(dlong))
    ndat[,timecol]=time
    ndat[,idcol]=idi
    ndat[match(odat[,timecol],ndat[,timecol]),]=odat
    out<-rbind(out,ndat)
  }
  colnames(out) <- colnames(dlong)
  
  l1 <- sum(!is.na(dlong[,-which(colnames(dlong) %in% c(timecol,idcol))])) 
   l2 <-   sum(!is.na(out[,-which(colnames(out) %in% c(timecol,idcol))])) 

if(l1!=l2) warning(paste0(l1-l2,' cells of data removed due to time overlap, consider reducing timestep'))
    
  out[,TDpredNames][is.na(out[,TDpredNames])] <- 0 
  out[,TIpredNames][is.na(out[,TIpredNames])] <- 0 
  
  

  return(out)
}
  
