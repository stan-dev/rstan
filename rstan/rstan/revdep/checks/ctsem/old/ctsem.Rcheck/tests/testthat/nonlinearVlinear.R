if( .Machine$sizeof.pointer != 4){
library(ctsem)
library(testthat)
set.seed(1)

context("nonlinearcheck") 

test_that("simplenonlinearcheck", { 
sunspots<-sunspot.year
 sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
 id <- 1
 time <- 1749:1924
datalong <- cbind(id, time, sunspots)

#setup model
 ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1, n.TDpred = 1,
  manifestNames='sunspots', 
  latentNames=c('ss_level', 'ss_velocity'),
   LAMBDA=matrix(c( 1, 'ma1'), nrow=1, ncol=2),
   DRIFT=matrix(c(0, 'a21', 1, 'a22'), nrow=2, ncol=2),
   TDPREDEFFECT=matrix(c('tdeffect',0),2),
   MANIFESTMEANS=matrix(c('mm'), nrow=1, ncol=1),
   MANIFESTVAR=diag(0,1),
   T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
   DIFFUSION=matrix(c(0, 0, 0, 'diff'), ncol=2, nrow=2))
 
 ssmodel$pars$indvarying<-FALSE #Because single subject
 ssmodel$pars$offset[14]<- 44 #Because not mean centered
 ssmodel$pars[4,c('transform','offset')]<- c(1,0) #To avoid multi modality 

 #td preds for testing only -- no real effect
 TD1 <- 0
 datalong <- cbind(datalong,TD1)
 datalong[seq(10,150,10),'TD1'] = 1
 
ssfitnl <- ctStanFit(datalong, ssmodel, iter=300, chains=2,optimize=T,verbose=0,maxtimestep = 2,
  nlcontrol=list(nldynamics=TRUE,ukffull=1),optimcontrol = list(finishsamples=50),nopriors=TRUE,deoptim=FALSE)
ssfitl <- ctStanFit(datalong, ssmodel, iter=300, chains=2,optimize=T,verbose=0,maxtimestep = 2,
  nlcontrol=list(nldynamics=FALSE),optimcontrol = list(finishsamples=50,deoptim=FALSE),nopriors=TRUE)

ssfitnlm <- ctStanFit(datalong, ssmodel, iter=300, chains=2,optimize=T,verbose=0,maxtimestep = 2,
  nlcontrol=list(nldynamics=TRUE,ukffull=1,nlmeasurement=TRUE),optimcontrol = list(finishsamples=50),nopriors=TRUE,deoptim=FALSE)

#output
snl=summary(ssfitnl)
snlm=summary(ssfitnlm)
sl=summary(ssfitl)

# expect_equal(snl$popmeans[,'mean'], sl$popmeans[,'mean'])
expect_equal(round(ssfitnl$stanfit$transformedpars_old[grep('pop_',rownames(ssfitnl$stanfit$transformedpars_old)),'mean'],3),
round(ssfitl$stanfit$transformedpars_old[grep('pop_',rownames(ssfitnl$stanfit$transformedpars_old)),'mean'],3),tol=.01)

expect_equal(round(ssfitnlm$stanfit$transformedpars_old[grep('pop_',rownames(ssfitnl$stanfit$transformedpars_old)),'mean'],3),
round(ssfitl$stanfit$transformedpars_old[grep('pop_',rownames(ssfitnl$stanfit$transformedpars_old)),'mean'],3),tol=.01)
})
}
