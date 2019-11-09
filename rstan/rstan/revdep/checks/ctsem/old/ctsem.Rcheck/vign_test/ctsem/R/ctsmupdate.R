ctsmupdate<-function(usecurrentwd=FALSE){
  
 sunspots<-datasets::sunspot.year
 sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
 id <- 1
 time <- 1749:1924
datalong <- cbind(id, time, sunspots)

#setup model
 model <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
  manifestNames='sunspots', 
  latentNames=c('ss_level', 'ss_velocity'),
   LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
   DRIFT=matrix(c(-.0001, 'a21', 1, 'a22'), nrow=2, ncol=2),
   MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
   CINT=matrix(c(0, 0), nrow=2, ncol=1),
   MANIFESTVAR=diag(.001,1),
   T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
   DIFFUSION=matrix(c(.0001, 0, 0, "diffusion"), ncol=2, nrow=2))

#fit
smgen <- ctStanFit(datalong, model,fit=FALSE,gendata=TRUE)$stanmodeltext
stanc(model_code = smgen,verbose = TRUE)

sm <- ctStanFit(datalong, model,fit=FALSE,gendata=FALSE)$stanmodeltext
stanc(model_code = sm,verbose = TRUE)

# model$w32 <- TRUE
# smgen32 <- ctStanFit(datalong, model,fit=FALSE,gendata=TRUE)$stanmodeltext
# stanc(model_code = smgen32,verbose = TRUE)
# 
# sm32 <- ctStanFit(datalong, model,fit=FALSE,gendata=FALSE)$stanmodeltext
# stanc(model_code = sm32,verbose = TRUE)

message(paste0('Update files? T / F?'))
continue <- readline()
if(continue){
  pathbase <- ifelse(usecurrentwd, paste0(getwd(),'/src/'),'~/../Dropbox/MPIB/CT-SEM/ctsem/src/')
  for(wi in 2){
    stan_files<-ifelse(wi==1,'stan_files32','stan_files')
  file.rename(paste0(pathbase,stan_files,'/ctsm.stan'), paste0('./src/',stan_files,'/ctsm.bak'))
  file.rename(paste0(pathbase,stan_files,'/ctsmgen.stan'), paste0('./src/',stan_files,'/ctsmgen.bak'))
sink(file=paste0(pathbase,stan_files,'/ctsm.stan'))
# if(wi==1) cat(sm32) else 
cat(sm)
sink()
sink(file=paste0(pathbase,stan_files,'/ctsmgen.stan'))
# if(wi==1) cat(smgen32) else 
cat(smgen)
sink()
}
# 
# message('All ok? finish this...')
# compile <- readline()
# if(compile) eval(parse(text=paste0('eval(devtools::install(local=FALSE),envir = globalenv())')))
}
}
