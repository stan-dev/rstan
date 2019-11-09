require(ctsem)
require(testthat)

context("kalmanVram")

test_that("time calc", {
  set.seed(4)
nsubjects=20
  #2 latent 2 td preds
TDPREDEFFECT=matrix(c(-20,1,-1,20),2,2)
  
gm=ctModel(Tpoints=10,n.latent=2,n.manifest=3,LAMBDA=matrix(c(1,0,0,0,1,.5),ncol=2),
  DRIFT=matrix(c(-.3,0,0,-.05),2,2),DIFFUSION=diag(2,2),
  # TRAITVAR=diag(2),
  CINT=matrix(c(5,3),ncol=1),
  TDPREDEFFECT=TDPREDEFFECT,
  TDPREDVAR=diag(.002,10*2),
  MANIFESTVAR=diag(.3,3),
  T0VAR=diag(2),
  n.TDpred=2,
  TDPREDMEANS=matrix(rep(c(1,rep(0,4),0,rep(0,4)),2),ncol=1))

gd=ctGenerate(gm,nsubjects,burnin=50)

m=ctModel(Tpoints=10,n.latent=2,n.manifest=3,LAMBDA=matrix(c(1,0,0,0,1,.5),ncol=2),
  DRIFT=matrix(c(-.3,0,0,-.05),2,2),
  # MANIFESTVAR=diag(.3,3),
  # DIFFUSION=diag(2,2),
  TDPREDEFFECT=matrix(c('td1',1,-1,20),2,2),
  # TRAITVAR='auto',
  CINT=matrix(c(5,3),ncol=1),
  MANIFESTMEANS=matrix(c(0,0,0),ncol=1),
  n.TDpred=2)

f1=ctRefineTo(dat = gd,m,retryattempts = 3,objective='Kalman',stationary=c('T0MEANS','T0VAR'),carefulFit=TRUE)
# f1$mxobj=mxRun(f1$mxobj)
f2=ctRefineTo(dat = gd,m,retryattempts = 3,objective='mxRAM',stationary=c('T0MEANS','T0VAR'),carefulFit=TRUE)

expect_equal(f1$mxobj$output$estimate[-1],f2$mxobj$output$estimate[-1],tolerance=.001)

ctPostPredict(f1,timestep=.1,n.subjects=200)

})
