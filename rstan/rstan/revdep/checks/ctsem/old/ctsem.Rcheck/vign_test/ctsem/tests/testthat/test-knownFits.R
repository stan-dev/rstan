library(ctsem)
library(testthat)

context("knownFits")



#anomauth
test_that("anomauth", {

  data(AnomAuth)
  AnomAuthmodel<-ctModel(LAMBDA=matrix(c(1, 0, 0, 1), nrow=2, ncol=2),  
    n.latent=2,n.manifest=2, 
    MANIFESTVAR=diag(0,2),
    Tpoints=5)

  ll=Inf
  counter=0
  while(counter < 10 && (ll >  23415.99  || ll < 23415.0)){
    counter=counter+1
  AnomAuthfit<-ctFit(AnomAuth, AnomAuthmodel)
  ll=AnomAuthfit$mxobj$output$Minus2LogLikelihood
  }
  
  expect_equal(23415.929,ll)
  
 if( .Machine$sizeof.pointer != 4){
   sm <- ctStanModel(AnomAuthmodel)
  sm$pars$indvarying<- FALSE
  sf=ctStanFit(ctDeintervalise(ctWideToLong(AnomAuth,Tpoints = AnomAuthmodel$Tpoints,n.manifest = 2)),
    ctstanmodel = sm, optimize=TRUE,verbose=0,savescores = FALSE,nopriors=TRUE,optimcontrol=list(finishsamples=10))
  expect_equal(23415.929,-2*sf$stanfit$optimfit$value,tolerance=.01)
 }

})

#anomauth with trait asymptotic vs standard param comparisons
test_that("time calc", {
  
  data(AnomAuth)
  AnomAuthmodel<-ctModel(LAMBDA=matrix(c(1, 0, 0, 1), nrow=2, ncol=2),  
    n.latent=2,n.manifest=2, 
    MANIFESTVAR=diag(0,2),
    TRAITVAR='auto',
    CINT=matrix(0,nrow=2),
    T0MEANS=matrix(0,nrow=2),
    MANIFESTMEANS=matrix(c('m1','m2'),nrow=2),
    timeVarying='MANIFESTMEANS',
    Tpoints=5)
  
  AnomAuthfit1<-ctRefineTo(AnomAuth, AnomAuthmodel,asymptotes=FALSE, verbose=0,retryattempts=1)
  AnomAuthfit2<-ctRefineTo(AnomAuth, AnomAuthmodel,asymptotes=TRUE, verbose=0,retryattempts=1)
  
  
  expect_equal(AnomAuthfit2$mxobj$output$Minus2LogLikelihood,AnomAuthfit1$mxobj$output$Minus2LogLikelihood)
  
  summ1<-summary(AnomAuthfit1,verbose=TRUE)
  summ2<-summary(AnomAuthfit2,verbose=TRUE)
  
  expect_equal(summ1$ctparameters[,1:2],summ2$ctparameters[,1:2],tolerance = .001)
  
})


test_that("oscillator", {
data("Oscillating")

inits <- c(-39.5, -.5, .1, 1, 0, 1, 0.05, .9)
names(inits) <- c("crosseffect","autoeffect", "diffusion",
  "T0var11", "T0var21", "T0var22","m1", "m2")

oscillatingm <- ctModel(n.latent = 2, n.manifest = 1, Tpoints = 11, 
  MANIFESTVAR = matrix(c(0), nrow = 1, ncol = 1), 
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  T0MEANS = matrix(c('m1', 'm2'), nrow = 2, ncol = 1), 
  T0VAR = matrix(c("T0var11", "T0var21", 0, "T0var22"), nrow = 2, ncol = 2),
  DRIFT = matrix(c(0, "crosseffect", 1, "autoeffect"), nrow = 2, ncol = 2), 
  CINT = matrix(0, ncol = 1, nrow = 2),
  DIFFUSION = matrix(c(0, 0, 0, "diffusion"), nrow = 2, ncol = 2))#,
  # startValues = inits)

oscillatingf <- ctFit(Oscillating, oscillatingm, carefulFit = FALSE,retryattempts = 3)

expect_equal(-3461.936,oscillatingf$mxobj$output$Minus2LogLikelihood,tolerance=.001)

if( .Machine$sizeof.pointer != 4){
 sm <- ctStanModel(oscillatingm)
  sm$pars$indvarying<- FALSE
  sf=ctStanFit(ctDeintervalise(ctWideToLong(Oscillating,Tpoints = oscillatingm$Tpoints,n.manifest = 1)),
    ctstanmodel = sm, optimize=TRUE,verbose=0,savescores = FALSE,nopriors=TRUE,optimcontrol=list(finishsamples=10))
  expect_equal(-3461.936,-2*sf$stanfit$optimfit$value,tolerance=.01)
  
}

})
