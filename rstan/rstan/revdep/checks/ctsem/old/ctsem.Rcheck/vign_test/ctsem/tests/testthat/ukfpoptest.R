if(Sys.getenv("NOT_CRAN")==TRUE & .Machine$sizeof.pointer != 4){

library(ctsem)
library(testthat)
set.seed(1)

context("ukfpopcheck") 

test_that("ukfpopcheck1", {
Tpoints<-20
n.latent=1
n.manifest=1
nsubjects=300
burnin=30
dtmat = matrix(exp(rnorm(burnin+Tpoints-1,-.3,.5)),1)
par1<-rnorm(nsubjects,-.3,.8)
for(i in 1:nsubjects){
gm<-ctModel(type='omx',n.latent=1,n.manifest=n.manifest,Tpoints=Tpoints,LAMBDA=matrix(rep(1,n.manifest),ncol=1),
  DRIFT=diag(-.4,1),
  CINT=diag(par1[i],1),
  T0VAR=diag(.1,1),
  # MANIFESTMEANS=diag(par1[i],1), #matrix(mmeans,nrow=n.manifest),
  MANIFESTVAR=t(chol(diag(.0001,n.manifest))),
  DIFFUSION=t(chol(diag(3,1))))
if(i==1) cd<-ctGenerate(gm,n.subjects=1,burnin=burnin,wide=FALSE,dtmat = dtmat) else {
  newdat <- ctGenerate(gm,n.subjects=1,burnin=burnin,wide=FALSE,dtmat = dtmat)
  newdat[,'id'] <- i
  cd<-rbind(cd,newdat)
}
}

cd[,gm$manifestNames]<-(cd[,gm$manifestNames]) + rnorm(length(cd[,gm$manifestNames]),0, .9^2)


m1<-ctModel(type='omx',n.latent=2,n.manifest=1,Tpoints=Tpoints,
  LAMBDA=matrix(c(1,1),nrow=n.manifest,ncol=2),
  MANIFESTMEANS=matrix(0,nrow=n.manifest),
  # CINT=matrix(c('cint1',0),2,1),
  # T0MEANS=matrix(c('t0m1',0),2),
  # T0VAR=matrix(c('t0var11',0, 0,'t0var22'),2,2),
  DIFFUSION=matrix(c('diff11',0,0,1e-5),2,2),
  DRIFT=matrix(c('dr11',0,1,-1e-5),2,2)
)
# m1$T0VAR[2,1]=0

#model for ukf ctsem
m2<-ctModel(type='omx',n.latent=1,n.manifest=n.manifest,Tpoints=Tpoints,
  # T0MEANS=matrix(c(0),1),
  # T0VAR=matrix(c(1.5),1),
  MANIFESTMEANS=matrix(0),CINT=matrix('cint'),
  # DIFFUSION=matrix(c(1.4),1),
  # DRIFT=matrix(c(-.4),1),
  TRAITVAR='auto',
    LAMBDA=matrix(1),
  )



# #original ctsem
cfit1 <- ctRefineTo(dat = cd,dataform = 'long',ctmodelobj = m1,retryattempts = 1)
cfit1$mxobj$DRIFT$values
cfit1$mxobj$DIFFUSION$result
cfit1$mxobj$T0VAR$result

cfit2 <- ctRefineTo(dat = cd,dataform = 'long',ctmodelobj = m2,retryattempts = 0,stationary='',carefulFit=T)
cfit2$mxobj$DRIFT$values
cfit2$mxobj$DIFFUSION$result
cfit2$mxobj$T0VAR$result

summary(cfit1)$ctparameters
summary(cfit2)$ctparameters

#bayesian / ukf ctsem
sm1 <- ctStanModel(m1)
sm1$pars$indvarying <- FALSE
# sm1$pars$indvarying[!sm1$pars$matrix %in% c('MANIFESTMEANS')] <- FALSE

# sink(file='../sf1.txt')
sf1 <- ctStanFit(cd,sm1,iter=200,chains=4,cores=4,
  optimize=TRUE,verbose=1,nopriors = TRUE,
  optimcontrol=list(deoptim = FALSE,isloops=0,isloopsize=50,finishsamples=50),
  derrind=1,
  nlcontrol=list(nldynamics=FALSE))
# sink()
summary(sf1)$popmeans
s1=sf1$stanfit$transformedpars_old[grep('pop_DRIFT',rownames(sf1$stanfit$transformedpars_old)),'mean'][1]


m2$T0VAR[,]=0
# m2$T0MEANS[,]=2.588
sm2 <- ctStanModel(m2)
sm2$pars$indvarying[!(sm2$pars$matrix %in% c('CINT','T0MEANS'))] <- FALSE

# sink(file='../sinkout.txt')
sf2 <- ctStanFit(cd,sm2,iter=200,chains=4,cores=4,
  optimize=TRUE,verbose=1,nopriors = TRUE,intoverpop = TRUE,
  optimcontrol=list(deoptim = FALSE,isloops=0,isloopsize=50,finishsamples=50),
  derrind=1)
# sink()

summary(sf2)$popmeans
s2=sf2$stanfit$transformedpars_old[grep('pop_DRIFT',rownames(sf2$stanfit$transformedpars_old)),'mean']

expect_equivalent(s1,s2,cfit2$mxobj$DRIFT$values,cfit1$mxobj$DRIFT$values,tol=1e-3)
})
}
