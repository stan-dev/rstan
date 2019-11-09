databasedTIpreds <- function(dlong,manifestNames, TIpredNames = NA, order=2){

  dat <- sapply(c('sd','mean'), function(f){
    func <- eval(parse(text=f))
    sapply(manifestNames, function(vari) {
      y<-cbind(sapply(unique(dlong[,'id']), function(idi){
        x<-func(dlong[dlong[,'id'] %in% idi,vari],na.rm=TRUE)
        names(x) <- idi
        return(x)
      },simplify = TRUE))
      colnames(y) <- paste0(vari,'_',f)
      return(y)
    },simplify = TRUE)
  },simplify = "array")


  dat2<-matrix(dat,nrow=dim(dat)[1])
  colnames(dat2) <- paste0('z_',dimnames(dat)[[2]], '_',rep(dimnames(dat)[[3]],each=length(dimnames(dat)[[2]])))
  if(!is.na(TIpredNames[1])) {
    scaledTI <-dlong[match(unique(dlong[,'id']),dlong[,'id']) ,TIpredNames,drop=FALSE]
    colnames(scaledTI) <- paste0('z_',colnames(scaledTI))
    originalTI <- dlong[match(unique(dlong[,'id']),dlong[,'id']) ,TIpredNames,drop=FALSE]
    dat2 <- cbind(dat2,scaledTI)
  }
  
  if(order > 1){
    for(oi in 2:order){
    dat2o <- scale(dat2)^order #scale and center before powering to reduce colinearity
    colnames(dat2o) <- paste0(colnames(dat2),'^',order)
    dat2 <- cbind(dat2,dat2o)
    }
  }
    
  
  TIpredNames=colnames(dat2)
  id<-unique(dlong[,'id'])
  dat3 <- cbind(id,scale(dat2),originalTI)
  
  
  dat3<-merge(x = dlong[,!colnames(dlong) %in% TIpredNames], y=dat3)
  return(list(dat=dat3, TIpredNames=TIpredNames))
}
  
  
