if(Sys.getenv("NOT_CRAN")==TRUE & .Machine$sizeof.pointer != 4){

library(ctsem)
library(testthat)
set.seed(1)

context("tipredcheck")

test_that("simpleTIpredcheck", {
Tpoints=100
n.manifest=1
n.TDpred=0
n.TIpred=1
n.latent=1
n.subjects=10
TI1 <- rnorm(n.subjects)
gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
n.TDpred=n.TDpred,n.manifest=n.manifest,
  MANIFESTVAR=diag(0.5,1),
  LAMBDA=diag(1,1),
  DRIFT=matrix(c(-.3),nrow=1),
  DIFFUSION=matrix(c(2),1),
  T0VAR=diag(10,1))

for(i in 1:n.subjects){
  gm$CINT[1,1] <- TI1[i]*5
ndat<-ctGenerate(gm,n.subjects=1,burnin=30,wide=FALSE,logdtsd=.4)
ndat <- cbind(ndat,TI1[i])
ndat[,1] <- i
if(i>1) tdat <- rbind(tdat,ndat) else tdat <- ndat
}
colnames(tdat)[4] <- 'TI1'

tdat[2,'Y1'] <- NA
tdat[tdat[,'id']==2,'TI1'] <- NA

checkm<-ctModel(type='stanct',Tpoints=Tpoints,
MANIFESTVAR=diag(0.5,1),
  DRIFT=matrix(c(-.3),nrow=1),
  DIFFUSION=matrix(c(2),1),
  n.latent=n.latent,n.TDpred=n.TDpred,n.TIpred=n.TIpred,
  MANIFESTMEANS=matrix(0,nrow=n.manifest),
  CINT=matrix(c('cint1'),ncol=1),
  n.manifest=n.manifest,LAMBDA=diag(1))

 checkm$pars$indvarying <- FALSE

 checkm$pars[c(-1,-7) ,c('TI1_effect')] <- FALSE
 checkm$pars$meanscale[c(1,7)] <- 1

tfit<-ctStanFit(tdat,checkm,iter=400, warmup=380,chains=2,optimize=T,isloops=3,nopriors=F,verbose=0)
s=summary(tfit)

expect_equal(s$tipreds[2,'mean'],5,tolerance=.1)

tfit<-ctStanFit(tdat,checkm,iter=400, warmup=380,chains=2,optimize=T,isloops=3,nopriors=F,verbose=0,
  nlcontrol=list(nldynamics=TRUE,nlmeasurement=TRUE))
s=summary(tfit)

expect_equal(s$tipreds[2,'mean'],5,tolerance=.1)

tfit<-ctStanFit(tdat,checkm,iter=400,chains=2,control=list(adapt_delta=.8,max_treedepth=6))
s=summary(tfit)

expect_equal(s$tipreds[2,'mean'],5,tolerance=.1)
})
}
