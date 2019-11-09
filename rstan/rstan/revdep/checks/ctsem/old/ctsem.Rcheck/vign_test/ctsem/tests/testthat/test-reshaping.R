library(ctsem)
library(testthat)

context("reshaping")

test_that("time calc", {

  Tpoints=3
  n.latent=2
  n.manifest=4
  n.TDpred=2
  n.TIpred=3
  
  generatingModel<-ctModel(Tpoints=Tpoints,n.latent=n.latent,
    n.TDpred=n.TDpred,n.TIpred=n.TIpred,n.manifest=n.manifest,
    LAMBDA=matrix(c(1,.4,.8,0,0,0,0,1),nrow=n.manifest,ncol=n.latent),
    MANIFESTVAR=diag(c(1),n.manifest),
    TRAITVAR=matrix(c(2.3,-1.1,0,1.8) ,n.latent,n.latent),
    MANIFESTMEANS=matrix(c(0,0,0,0),n.manifest,1),
    DRIFT=matrix(c(-.23,.1,.0,-.4),n.latent),
    DIFFUSION=matrix(c(8.3,-5.1,0,8.4),n.latent,n.latent),
    CINT=matrix(c(0,.4),n.latent,1),
    TRAITTDPREDCOV = matrix(c(.6,-.3,.4,.4),nrow=n.latent,ncol=n.TDpred*(Tpoints)),
    TDPREDEFFECT=matrix(c(1.2,-.4, 0,.3),nrow=n.latent,ncol=n.TDpred),
    TIPREDEFFECT=matrix(c(.32,.08,.4,-.6,-.3,0),nrow=n.latent,ncol=n.TIpred),
    # TDPREDMEANS=matrix(c(0,0,1,rep(0, (Tpoints-1-3)*n.TDpred)),nrow=n.TDpred*(Tpoints-1)),
    TIPREDMEANS=matrix(0,nrow=n.TIpred),
    TDPREDVAR=diag(1,n.TDpred*(Tpoints)),
    TIPREDVAR=diag(.5,n.TIpred),
    T0MEANS=matrix(0,ncol=1,nrow=n.latent))
  
  data<-ctGenerate(generatingModel,n.subjects=20,burnin=500)
  data[1:(ceiling(nrow(data/2))),'dT1']<-2
  data[1:(ceiling(nrow(data/2))),'dT2']<-3

  
  manifestNames<-paste0('manifestV',1:n.manifest)
  latentNames<-paste0('latentV',1:n.latent)
  TDpredNames<-paste0('TDpredV',1:n.TDpred)
  TIpredNames<-paste0('TIpredV',1:n.TIpred)
  
  colnames(data)<-ctWideNames(n.manifest=n.manifest, n.TDpred = n.TDpred,
    Tpoints=Tpoints, manifestNames=manifestNames, TDpredNames=TDpredNames, TIpredNames=TIpredNames, n.TIpred=n.TIpred)
  
  testlong<-ctWideToLong(data,n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred,
    manifestNames=manifestNames,TDpredNames=TDpredNames,TIpredNames=TIpredNames,Tpoints=Tpoints)
  
  testlong<-ctDeintervalise(testlong)
  
  testwide<-ctLongToWide(testlong,id='id',time='time',
    manifestNames=manifestNames,TDpredNames=TDpredNames,TIpredNames=TIpredNames)
  
  testwide<-ctIntervalise(testwide,n.manifest=n.manifest,n.TDpred=n.TDpred,n.TIpred=n.TIpred, Tpoints=Tpoints,
    manifestNames=manifestNames,TDpredNames=TDpredNames,TIpredNames=TIpredNames)
  
  identical(testwide,data)
  

  expect_identical(testwide,data)
})



