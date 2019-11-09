priorchecker <- function(sf,pars=c('rawpopmeans','rawpopsdbase','tipredeffectparams'),digits=2){
  e=extract(sf)
  funcs <- c(base::mean,stats::sd)
  pars=unlist(lapply(pars,function(x) if(!is.null(dim(e[[x]]))) x))
  out=round(do.call(cbind,lapply(funcs, function(fn) do.call(c,
    lapply(pars, function(obji) apply(e[[obji]],2,fn,na.rm=TRUE)) ))),digits)
  rownames(out)=do.call(c,c(lapply(pars, function(obji) paste0(obji,'_',1:ncol(e[[obji]])))))
  out=data.frame(out,do.call(c,c(lapply(pars, function(obji) 1:ncol(e[[obji]])))))
  out=data.frame(out,do.call(c,c(lapply(pars, function(obji) rep(obji, ncol(e[[obji]]))))),stringsAsFactors = FALSE)
  colnames(out)=c('mean','sd', 'param', 'object')
  rownames(out) = getparnamesfromraw(priorcheck=out,sf=sf)
  return(out)
}

getparnamesfromraw <- function(priorcheck, sf){
  newnames=rownames(priorcheck)
  for(ni in 1:nrow(priorcheck)){
    if(priorcheck$object[ni] %in% 'rawpopmeans'){
      newnames[ni]=paste0('rawpop_',sf$setup$popsetup$parname[sf$setup$popsetup$param %in% priorcheck$param[ni]][1])
    }
    if(priorcheck$object[ni] %in% 'tipredeffectparams'){
      newnames[ni]=paste0('rawtipredeffect_',paste0(
        which(sf$standata$TIPREDEFFECTsetup == priorcheck$param[ni],arr.ind = TRUE),collapse='_'))
    }
  }
  return(newnames)
}

priorcheckreport <- function(sf, meanlim = 2, sdlim= .2,digits=2){
  p=priorchecker(sf)
  ps=sf$setup$popsetup
  p=p[abs(p$mean) > meanlim | p$sd > sdlim,]
  out<-list(priorcheck_note='The following posteriors exceeded arbitrary limits re normal(0,1) -- priors / transforms are likely somewhat informative. Not necessarily a problem.')
  out$priorcheck=p[,c('mean','sd')]
  
  # if(any(p$object %in% 'rawpopsdbase')){
  #   e=apply(extract(sf,pars='rawpopsdprops')$rawpopsdprops,2,mean,na.rm=TRUE)
  #   names(e) = ps$parname[match(x = 1:length(e),ps$param)]
  #   e=e[e> 1/length(e) | e==max(e)]
  #   out$priorcheck_sd_note = 'Population posterior variance exceeded check limits. Not necessarily a problem, but these parameters contribute most variance: '
  #   out$priorcheck_sd = round(e,digits)
  # }
  return(out)
}

