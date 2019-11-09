ctdataupdate<-function(forcerecompile=FALSE){
  Tpoints=20
  n.manifest=2
  n.TDpred=1
  n.TIpred=3
  n.latent=2
  n.subjects=5
  gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
    n.TDpred=n.TDpred,n.TIpred=n.TIpred,n.manifest=n.manifest,
    MANIFESTVAR=diag(0.5,2),
    TIPREDEFFECT=matrix(c(.5,0,0,-.7,0,2),nrow=2),
    TIPREDVAR=matrix(c(1,-.2,0, 0,1,0, 0,0,.5),nrow=3),
    TDPREDVAR=matrix(0,nrow=n.TDpred*(Tpoints),ncol=n.TDpred*(Tpoints)),
    TDPREDMEANS=matrix(round(exp(rnorm(n.TDpred*(Tpoints),-1.9,1)),0),
      nrow=n.TDpred*(Tpoints)),
    TDPREDEFFECT = matrix(c(1,-1),ncol=1),
    LAMBDA=diag(1,2),
    DRIFT=matrix(c(-.3,.2,0,-.2),nrow=2),
    DIFFUSION=matrix(c(2,1,0,2),2),
    CINT=matrix(c(0,0),nrow=2),
    T0MEANS=matrix(0,ncol=1,nrow=2),
    T0VAR=diag(100,2))
  
  checkm<-ctModel(type='stanct',
    n.latent=n.latent,n.TDpred=n.TDpred,n.TIpred=n.TIpred,
    MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
    MANIFESTMEANS=matrix(0,nrow=n.manifest),
    CINT=matrix(c('cint1','cint2'),ncol=1),
    n.manifest=n.manifest,LAMBDA=diag(2))
  
  checkm$pars$indvarying[-c(7,13)] <- FALSE
  checkm$pars$sdscale <- .1
  
  checkm$pars[c(-1,-2, -21,-22) ,c('TI1_effect','TI2_effect','TI3_effect')] <- FALSE
  
  ctstantestfit<-ctStanFit(ctstantestdat,checkm,iter=500, warmup=460,thin=2,chains=2,
    forcerecompile=forcerecompile,
    control=list(max_treedepth=8,adapt_delta=.8),save_warmup=FALSE)
  ctstantestfit <- ctStanGenerateData(ctstantestfit)
  summary(ctstantestfit)
  message(paste0('Updating from ',(getwd()),', continue T / F?'))
continue <- readline()
if(continue){
  save(ctstantestfit,file='.\\data\\ctstantestfit.rda')
  paths <- sort(Sys.glob(c("data/*.rda", "data/*.RData")))
  tools::resaveRdaFiles(paths)
}
}
