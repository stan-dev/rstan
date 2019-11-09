#' Get Kalman filter estimates from a ctStanFit object
#'
#' @param fit fit object from \code{\link{ctStanFit}}.
#'
#' @return list containing Kalman filter elements, each element in array of
#' iterations, data row, variables. llrow is the log likelihood for each row of data.
#' @export
#'
#' @examples 
#' ctStanKalman(ctstantestfit)
#' 
ctStanKalman <- function(fit){
  k=extract(fit)$kalman
  k[k==99999] <- NA #for missingness
  nlatent <- fit$standata$nlatentpop
  nmanifest <- fit$standata$nmanifest
  dimnames(k) = list(iter=1:dim(k)[1],drow=1:dim(k)[2],
    kalman=paste0(c(rep('lln',nmanifest),
      rep('llscale',nmanifest),rep('err',nmanifest),rep('ypred',nmanifest),rep('etaprior',nlatent),rep('etaupd',nlatent)),
      c(1:nmanifest,1:nmanifest,1:nmanifest,1:nmanifest,1:nlatent,1:nlatent)))
  
  
  
  obj <- c('lln','llscale','err','ypred','etaprior','etaupd')
  
  
  lln=k[,,1:nmanifest,drop=FALSE]
  llscale=k[,,(nmanifest*1+1):(nmanifest*1+nmanifest),drop=FALSE]
  err=k[,,(nmanifest*2+1):(nmanifest*2+nmanifest),drop=FALSE]
  ypred=k[,,(nmanifest*3+1):(nmanifest*3+nmanifest),drop=FALSE]
  etaprior=k[,,(nmanifest*4+1):(nmanifest*4+nlatent),drop=FALSE]
  etaupd=k[,,(nmanifest*4+nlatent+1):(nmanifest*4+nlatent*2),drop=FALSE]
  
  llvec = apply(lln,1:2,function(x) {
    sum(dnorm(x[!is.na(x)],log = TRUE))
  })
  llrow = llvec - apply(llscale, 1:2, function(x) sum(x,na.rm=TRUE))

    return(list(lln=lln,llscale=llscale,err=err,ypred=ypred,etaprior=etaprior,etaupd=etaupd,llrow=llrow))
}



