#current disabled in ctsem
if(1==0){
ctcorplots <- function(dlong, vars, maxlag, splitvar=NA, splitpoints=NA,varsPerPlot=Inf,splitAlpha=.3){
  dlong <- data.table(dlong)
  for(lead in 0:maxlag){
    nm1 <- vars
    nm2 <- paste("lead",lead,nm1, sep=".")
    dlong[, (nm2) := shift(.SD, type='lead',n = lead), by = id, .SDcols=nm1]
  }
  
  covvars <- c()
  for(vi in 1:length(vars)){
    covvars <- c(covvars,vars[vi])
    for(li in 1:maxlag){
      covvars <- c(covvars,paste('lead',li,vars[vi],sep='.'))
    }}
  
  nlcor <- FALSE
  corm <- list()
  if(is.na(splitvar) || all(is.na(splitpoints))) splitpoints <- c(-Inf,Inf) else splitpoints <- c(-Inf,splitpoints,Inf)
  for(si in 2:(length(splitpoints))){
    if(!(is.na(splitvar) || all(is.na(splitpoints)))) selection = dlong[[splitvar]]>splitpoints[si-1] & dlong[[splitvar]] <= splitpoints[si] else selection = 1:nrow(dlong)
    if(!nlcor){
      corm[[si-1]]=cor(dlong[selection,covvars,with=FALSE],use='pairwise.complete.obs')
    } else {
      # library('infotheo')
      # corm[[si]]=mutinformation(dlong[dlong$fromgesis==0,covvars,with=FALSE])
    }
  }
  splitpoints=splitpoints[c(-1,-length(splitpoints))]
  mincor=min(0,sapply(corm,min,na.rm=TRUE))
  nvarsp=min(varsPerPlot,length(vars))
  for(v in vars){
    for(varseti in 1:ceiling(length(vars)/nvarsp)){
      smlvars <- vars[((varseti-1)*nvarsp+1 ) : min(length(vars),(varseti)*nvarsp)]
      for(si in 1:length(corm)){
        ms<-matrix(NA, maxlag+1,length(smlvars))
        for(li in 0:maxlag){
          for(ci in 1:length(smlvars)){
            varref = ci + ((varseti-1)*nvarsp)
            ms[li+1, ci] = corm[[si]][v,ifelse(li > 0, paste('lead',li,smlvars[ci], sep='.'),smlvars[ci])]
          }
        }
        colnames(ms) = smlvars
        colseq=0:(length(smlvars)-1) %% 5 +1
        varseq = order(ms[1,],decreasing = TRUE)
        
        colseqa=colseq
        ltyseq=1:length(smlvars)
        if(si > 1) {
          colseqa = sapply(colseq,function(x) adjustcolor(x, alpha.f = splitAlpha))
          ltyseq=1
        }
        if(si==1){
          plot.new()
          abline(v = 0, col = "grey80", lwd = 8000)
          par(new=TRUE)
        }
        matplot(seq(0,maxlag,1),ms[,varseq],type='b',main=v,col=colseqa,xlim=c(0,maxlag),ylim=c(mincor,1),lty = ltyseq,xlab='Lag',
          ylab=paste0(ifelse(nlcor,'MI','Corr.')),pch=1:length(smlvars),add = ifelse(si>1,TRUE,FALSE))
        grid(col = 'white')
        # matpoints(seq(0,to=(maxlag)/6, length.out=maxlag+1),mg[,varseq],col=colseqa,type='b',main=v,lty = 1:length(vars),lwd=2,pch=1:length(vars))
        if(si == length(corm)) legend('topright',legend = smlvars[varseq],text.col=colseq,lty=1:length(smlvars),bty='n',col=colseq,pch=1:length(smlvars))
      }
    }
  }
}



ctStanDataOut <- function(fit){
  dat <- data.table(fit$standata$time)
  setnames(dat, fit$ctstanmodel$timeName)
  dat[,(fit$ctstanmodel$subjectIDname) := data.table(fit$setup$idmap[fit$data$subject,1] ) ]
  dat[,(fit$ctstanmodel$manifestNames) := data.table(fit$data$Y)]
  if(fit$ctstanmodel$n.TDpred > 0) dat[,(fit$ctstanmodel$TDpredNames) := data.table(fit$standata$tdpreds)]
  if(fit$ctstanmodel$n.TIpred > 0) {
    tipreds = data.table(fit$setup$idmap[,1], fit$data$tipredsdata)
    setnames(tipreds, c(fit$ctstanmodel$subjectIDname,fit$ctstanmodel$TIpredNames ))
  }
  return(dat)
}

compareCor <- function(fit,maxlag,N,...){
  vars=fit$ctstanmodel$manifestNames
  newdat=data.table(ctStanDataOut(fit))
  dl2=data.table(newdat,new=0)
  gendat=ctStanGenerateData(fit,nsamples = N)$generate
  for(n in 1:N){
    ny=gendat$Y[,n,]
    newdat[,(vars) := data.table(ny)]
    newdat[,id:=paste0('n',N,'_',id)]
    dl2=rbind(dl2,data.table(newdat,new=n))
  }
  ctcorplots(dl2,vars = vars,maxlag = maxlag,splitvar = 'new',splitpoints = seq(1,N,1)-.5,...)
}

mixedcorplot <- function(cov1,cov2,corr=TRUE){
  if(corr) {
    bwvar <- cov2cor(cov1)
    # colnames(bwvar) <- vars
    # rownames(bwvar) <- vars
    bwvar[lower.tri(bwvar)] <- cov2cor(cov2)[lower.tri(bwvar)]
  } else {
    bwvar <- (cov1)
    bwvar[lower.tri(bwvar)] <- (cov2)[lower.tri(bwvar)]
  }
  # library(corrplot)
  corrplot(corr=bwvar,is.corr = corr)
}

} #end disabling
