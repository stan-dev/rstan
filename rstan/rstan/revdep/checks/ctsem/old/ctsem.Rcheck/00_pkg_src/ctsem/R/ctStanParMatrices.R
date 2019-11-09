#' Returns population system matrices from a ctStanFit object, and vector of values for free parameters.
#'
#' @param fit ctStanFit object.
#' @param parvalues vector of parameter values to assign to free parameters in the model
#' @param timeinterval time interval to use for discrete time (dt) matrix calculations.
#' @param sf stanfit object. Generally not necessary, but for repeated calls to this function, can speed things up.
#'
#' @return A list containing various matrices related to a continuous time dynamic model. 
#' Matrices with "dt" in front refers to discrete time, "asym" refers to asymptotic (time interval = infinity), 
#' and "cor" refers to correlations. 
#' @export
#'
#' @examples
#' \dontrun{
#' ctStanParMatrices(ctstantestfit,rnorm(17,0,.1))
#' }
ctStanParMatrices <- function(fit, parvalues, timeinterval=1, sf=NA){
  
  if(class(fit) !='ctStanFit') stop('not a ctStanFit object')
  model <- fit$ctstanmodel
  fit$standata$savescores <- 0L
  fit$standata$gendata <- 0L
  fit$standata$dokalman <- 0L
  # browser()
  if(length(parvalues)!=fit$data$nparams) stop('length of parvalues != number of free params (',fit$data$nparams,') in model!')
  if(suppressWarnings(is.na(sf))) sf <- stan_reinitsf(fit$stanmodel,data=fit$standata) #suppressOutput(sf <- suppressWarnings(sampling(,iter=1,control=list(max_treedepth=1),chains=1)))
  npars <- get_num_upars(sf)
  pars <- c(parvalues,rep(0,npars - fit$data$nparams))
  sfc <- constrain_pars(sf, pars)
  
  
   whichmatrices='all'
  
  if(whichmatrices[1] == 'all') {
    whichmatrices <- c(fit$setup$matrices$base,
      'asymDIFFUSION','dtDRIFT','dtDIFFUSION','dtCINT','DIFFUSIONcor','asymDIFFUSIONcor','T0VARcor','asymCINT')
  } else whichmatrices <- unique(c(whichmatrices, fit$setup$matrices$base)) #need base matrices for computations
  
  stanmats <- c(fit$setup$matrices$base,'asymDIFFUSION')[c(fit$setup$matrices$base,'asymDIFFUSION') %in% whichmatrices]
  
  out <- list()
  for(m in c(stanmats)){
    # assign(m,ctCollapse(sfc[[paste0('pop_',m)]],1,mean)) 
    out[[m]] <- sfc[[paste0('pop_',m)]]
  }
  
  #dimension naming (latent row object, manifest column object, etc
  for(lro in c('DRIFT','DIFFUSION','CINT','T0VAR','T0MEANS','asymDIFFUSION',if('TDPREDEFFECT' %in% model$pars$matrix) 'TDPREDEFFECT')){
    if(lro %in% whichmatrices) rownames(out[[lro]]) <- model$latentNames
  }
  for(lco in c('DRIFT','DIFFUSION','T0VAR','asymDIFFUSION','LAMBDA')){
    if(lco %in% whichmatrices) colnames(out[[lco]]) <- model$latentNames
  }
  for(mro in c('LAMBDA','MANIFESTVAR','MANIFESTMEANS')){
    if(mro %in% whichmatrices) rownames(out[[mro]]) <- model$manifestNames
  }
  for(mco in c('MANIFESTVAR')){
    if(mco %in% whichmatrices) colnames(out[[mco]]) <- model$manifestNames
  }
  
  if('TDPREDEFFECT' %in% model$pars$matrix)  colnames(out$TDPREDEFFECT) <- model$TDpredNames
  
  
  choltrue <- FALSE #!as.logical(fit$data$lineardynamics)
  
  # if(choltrue) DIFFUSION = msquare(DIFFUSION) #sdcovchol2cov(DIFFUSION,0)
  if('DIFFUSIONcor' %in% whichmatrices){
    out$DIFFUSIONcor = suppressWarnings(stats::cov2cor(out$DIFFUSION))
    out$DIFFUSIONcor[is.na(out$DIFFUSIONcor)] <- 0
  }
  
  
  # DRIFTHATCH<-DRIFT %x% diag(nrow(DRIFT)) + diag(nrow(DRIFT)) %x% DRIFT
  # asymDIFFUSION<-matrix(-solve(DRIFTHATCH, c(DIFFUSION)), nrow=nrow(DRIFT))
  if('asymDIFFUSIONcor' %in% whichmatrices){
    out$asymDIFFUSIONcor = suppressWarnings(stats::cov2cor(out$asymDIFFUSION))
    out$asymDIFFUSIONcor[is.na(out$asymDIFFUSIONcor)] <- 0
  }
  
  
  
  # dimnames(DRIFT)=list(ln,ln)
  # dimnames(DIFFUSION)=list(ln,ln)
  # dimnames(asymDIFFUSION)=list(ln,ln)
  # rownames(CINT)=ln
  # rownames(MANIFESTMEANS)=mn
  # dimnames(T0VAR)=list(ln,ln)
  # rownames(T0MEANS)=ln
  # dimnames(asymDIFFUSION)=list(ln,ln)
  # dimnames(LAMBDA)=list(mn,ln)
  # dimnames(MANIFESTVAR)=list(mn,mn)
  # out$MANIFESTVAR=MANIFESTVAR
  
  out$dtDRIFT=expm(out$DRIFT * timeinterval)
  if('dtDIFFUSION' %in% whichmatrices) out$dtDIFFUSION = out$asymDIFFUSION - (out$dtDRIFT %*% out$asymDIFFUSION %*% t(out$dtDRIFT ))
  if('dtDIFFUSIONcor' %in% whichmatrices) out$dtDIFFUSIONcor = cov2cor(out$dtDIFFUSION)
  if('dtCINT' %in% whichmatrices) out$dtCINT = (solve(out$DRIFT) %*%(out$dtDRIFT - diag(nrow(out$DRIFT))) %*% (out$CINT))
  if('asymCINT' %in% whichmatrices) out$asymCINT = -solve(out$DRIFT) %*% out$CINT
  
  if('asymDIFFUSIONcor' %in% whichmatrices) {
    out$T0VARcor = suppressWarnings(stats::cov2cor(out$T0VAR))
    out$T0VARcor[is.na(out$T0VARcor)] <- 0
  }
  
  
  # for(i in row(statspec)){
  #   if(statspec$matrix[i] =='T0VAR') {
  #     eval(parse(text=paste0(statspec$matrix[i], '[',statspec$row[i],' ,', statspec$col[i], '] <- ', 
  #     'asymDIFFUSION[',statspec$row[i],' ,', statspec$col[i], ']')))
  #     eval(parse(text=paste0(statspec$matrix[i], 'cor[',statspec$row[i],' ,', statspec$col[i], '] <- ', 
  #       'asymDIFFUSIONcor[',statspec$row[i],' ,', statspec$col[i], ']')))
  #   }
  #   if(statspec$matrix[i] =='T0MEANS') eval(parse(text=paste0(statspec$matrix[i], '[',statspec$row[i],' ,', statspec$col[i], '] <- ', 
  #     'asymCINT[',statspec$row[i],' ,', statspec$col[i], ']')))
  # }
  
  # T0VAR[upper.tri(T0VAR)] = t(T0VAR)[upper.tri(T0VAR)]
  # T0VARcor[upper.tri(T0VAR)] = t(T0VARcor)[upper.tri(T0VAR)]
  
  
  # out<-list(DRIFT=DRIFT,dtDRIFT=dtDRIFT, T0VAR=T0VAR, T0VARcor=T0VARcor, 
  #   DIFFUSION=DIFFUSION, DIFFUSIONcor=DIFFUSIONcor, dtDIFFUSION=dtDIFFUSION, dtDIFFUSIONcor=dtDIFFUSIONcor,
  #   asymDIFFUSION=asymDIFFUSION, asymDIFFUSIONcor=asymDIFFUSIONcor, 
  #   CINT=CINT,dtCINT=dtCINT, asymCINT=asymCINT, T0MEANS=T0MEANS,
  #   MANIFESTMEANS=MANIFESTMEANS, LAMBDA=LAMBDA)
  # 
  # 
  # if(choltrue) MANIFESTVAR=msquare(MANIFESTVAR)
  # 
  # 
  # 
  # if('TDPREDEFFECT' %in% model$pars$matrix) {
  #   dimnames(TDPREDEFFECT)=list(ln,tdn)
  #   out$TDPREDEFFECT<-TDPREDEFFECT
  # }
  
  return(out)
}
