
ctFitR<-function(datalong, ctmodel,...){
  
  carefulfit=FALSE
  regulariser=1
  
  # ctmodel <- list(ctmodel)
  kparsskeleton <- list(DRIFT=ctmodel$DRIFT,DIFFUSION=ctmodel$DIFFUSION,MANIFESTVAR=ctmodel$MANIFESTVAR,CINT=ctmodel$CINT,
    T0VAR=ctmodel$T0VAR,T0MEANS=ctmodel$T0MEANS,MANIFESTMEANS=ctmodel$MANIFESTMEANS,LAMBDA=ctmodel$LAMBDA)
  if(ctmodel$n.TDpred > 0) kparsskeleton$TDPREDEFFECT=ctmodel$TDPREDEFFECT
  if(!is.null(ctmodel$PARMEANS)) {
    kparsskeleton$PARMEANS <- ctmodel$PARMEANS
    kparsskeleton$PARVAR<- ctmodel$PARVAR
  }


  kpars <- unlist(utils::as.relistable(kparsskeleton))

  kparsfree <- which(suppressWarnings(is.na(as.numeric(kpars))))
  parnames <- unique(kpars[kparsfree])

  if(!is.null(ctmodel$calcs)) kparsskeleton$calcs <- ctmodel$calcs
  if(!is.null(ctmodel$manifesttype)) kparsskeleton$manifesttype <- ctmodel$manifesttype
  
  # update kpars with deterministic calcs
# browser()
# 
#   kparsskeleton$calcindices <- matrix(which(suppressWarnings(is.na(as.numeric(kpars))) & grepl('[',kpars,fixed=TRUE)),ncol=1)
#   kparsskeleton$calcs <- unique(kpars[kparsskeleton$calcindices])
#   kparsskeleton$calcmatrices <- unlist(lapply(kparsskeleton$calcs, function(x)
    
  kpars <- unlist(utils::as.relistable(kparsskeleton))

  npars <- length(parnames)
  parpositions <- unlist(lapply(kpars[kparsfree],function(x) which(parnames==x)))
  
  optimfunc <- function(pars,...){
    kpars[kparsfree] <- pars[parpositions]
    # PARMEANS <- as.numeric(kpars[kparsfree][names(kpars[kparsfree])=='PARMEANS'])
    # kpars[kparsskeleton$calcindices] <- unlist(lapply(kparsskeleton$calcs,function(x) eval(parse(text=x))))
    kpars <- utils::relist(kpars,kparsskeleton)
    
    
    if(runif(1) > .95)  print(matrix(c(parnames,pars),ncol=2))
    kpars[unlist(lapply(kpars,is.matrix))] <- lapply(kpars[unlist(lapply(kpars,is.matrix))],
      function(x) matrix(as.numeric(x),nrow=nrow(x),ncol=ncol(x)))
    kpars$DIFFUSION <- kpars$DIFFUSION %*% t(kpars$DIFFUSION)
    kpars$T0VAR <- kpars$T0VAR %*% t(kpars$T0VAR)
    kpars$MANIFESTVAR <- kpars$MANIFESTVAR %*% t(kpars$MANIFESTVAR)
    if(!is.null(kpars$PARVAR)) kpars$PARVAR <- kpars$PARVAR %*% t(kpars$PARVAR)
    kpars$DRIFT[c(row(kpars$DRIFT) == col(kpars$DRIFT)) & suppressWarnings(is.na(as.numeric(kparsskeleton$DRIFT)))] <- 
    -log(exp(-kpars$DRIFT)+1)[c(row(kpars$DRIFT) == col(kpars$DRIFT)) & suppressWarnings(is.na(as.numeric(kparsskeleton$DRIFT)))]
  

  lp=try(-Kalman(kpars=kpars,datalong=datalong,manifestNames=ctmodel$manifestNames,
    latentNames=ctmodel$latentNames,TDpredNames=ctmodel$TDpredNames,timecol='time',optimize=TRUE,...))
  if(is.nan(lp) | is.infinite(lp) | class(lp)=='try-error') lp <- -99999999999999
  if(carefulfit) lp <- lp +abs(-1 -(sum(diag(kpars$DRIFT))))*regulariser +abs(1 -(sum(diag(kpars$DIFFUSION))))*regulariser*.01
  return(lp)
  }
  
  fit <- optim(par= rnorm(npars,0,.5),fn=optimfunc,hessian=TRUE,
    method = c("L-BFGS-B"),control=list(trace=3,REPORT=1,
      ndeps=rep(1e-8,npars)),...)

  kpars[kparsfree] <- fit$par[parpositions]
  # PARMEANS <- as.numeric(kpars[kparsfree][names(kpars[kparsfree])=='PARMEANS'])
  # kpars[kparsskeleton$calcindices] <- unlist(lapply(kparsskeleton$calcs,function(x) eval(parse(text=x))))
  
  kpars <- utils::relist(kpars)

  kpars[unlist(lapply(kpars,is.matrix))] <- lapply(kpars[unlist(lapply(kpars,is.matrix))],
    function(x) matrix(as.numeric(x),nrow=nrow(x),ncol=ncol(x)))
  
  kpars$DIFFUSION <- kpars$DIFFUSION %*% t(kpars$DIFFUSION)
  kpars$T0VAR <- kpars$T0VAR %*% t(kpars$T0VAR)
  kpars$MANIFESTVAR <- kpars$MANIFESTVAR %*% t(kpars$MANIFESTVAR)

  kpars$DRIFT[c(row(kpars$DRIFT) == col(kpars$DRIFT)) & suppressWarnings(is.na(as.numeric(kparsskeleton$DRIFT)))] <- 
    -log(exp(-kpars$DRIFT)+1)[c(row(kpars$DRIFT) == col(kpars$DRIFT)) & suppressWarnings(is.na(as.numeric(kparsskeleton$DRIFT)))]
  
  # for(calci in kpars$calcs){
  #   eval(parse(text=calci))
  # }
  
  est <- unlist(kpars)[kparsfree]
  est <- est[c(names(est)[!duplicated(parpositions)])]
  stderrors <- try(sqrt(diag(solve(fit$hessian+diag(.000001,nrow(fit$hessian))))))
  if(length(stderrors)!=npars) stderrors <- rep(NA,npars)
  out <- matrix(c(est,stderrors),ncol=2)
  rownames(out) <- parnames
  colnames(out) <- c('Estimate', 'SE')
  return(out)
}
