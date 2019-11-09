if(.Machine$sizeof.pointer != 4){
library(ctsem)
library(testthat)

context("knownFits")



# #anomauth
# test_that("time calc", {
# 
#   data(AnomAuth)
#   AnomAuthmodel<-ctModel(type='omx',LAMBDA=matrix(c(1, 0, 0, 1), nrow=2, ncol=2),  
#     n.latent=2,n.manifest=2, 
#     MANIFESTVAR=diag(0,2),
#     Tpoints=5)
# 
#   AnomAuthfit<-ctFit(AnomAuth, AnomAuthmodel)
#   
#   dl<-ctWideToLong(datawide=AnomAuth,Tpoints=5,n.manifest=2)
#   dl<-ctDeintervalise(dl)
#   AnomAuthstanmodel<-ctModel(type='stanct',LAMBDA=matrix(c(1, 0, 0, 1), nrow=2, ncol=2),  
#     n.latent=2,n.manifest=2, 
#     MANIFESTVAR=diag(0.001,2),
#     Tpoints=5)
#   AnomAuthstanmodel$pars$indvarying<-FALSE
#   sfit<-ctStanFit(datalong=dl,ctstanmodel=AnomAuthstanmodel,
#     control=list(max_treedepth=6,adapt_delta=.8),
#     iter=500,chains=2,optimize=T,nopriors=T,verbose=1)
#   
#   # expect_equal(,AnomAuthfit$mxobj$DRIFT$values[1,1])
#   
# 
# })
# 
# #anomauth with trait asymptotic vs standard param comparisons
# test_that("time calc", {
#   
#   data(AnomAuth)
#   AnomAuthmodel<-ctModel(LAMBDA=matrix(c(1, 0, 0, 1), nrow=2, ncol=2),  
#     n.latent=2,n.manifest=2, 
#     MANIFESTVAR=diag(0,2),
#     TRAITVAR='auto',
#     CINT=matrix(0,nrow=2),
#     T0MEANS=matrix(0,nrow=2),
#     MANIFESTMEANS=matrix(c('m1','m2'),nrow=2),
#     timeVarying='MANIFESTMEANS',
#     Tpoints=5)
#   
#   AnomAuthfit1<-ctFit(AnomAuth, AnomAuthmodel,asymptotes=FALSE, objective='mxFIML',verbose=2,retryattempts=1)
#   AnomAuthfit2<-ctFit(AnomAuth, AnomAuthmodel,asymptotes=TRUE, verbose=2,retryattempts=1)
#   
#   
#   expect_equal(AnomAuthfit2$mxobj$output$Minus2LogLikelihood,AnomAuthfit1$mxobj$output$Minus2LogLikelihood)
#   
#   summ1<-summary(AnomAuthfit1,verbose=TRUE)
#   summ2<-summary(AnomAuthfit2,verbose=TRUE)
#   
#   expect_equal(summ1$ctparameters[,1:2],summ2$ctparameters[,1:2],tolerance = .001)
#   
# })


test_that("time calc", {
data("Oscillating")

# inits <- c(-39, -.5, 1, 1, .1, 1, 0, .9)
# names(inits) <- c("crosseffect","autoeffect", "diffusion",
#   "T0var11", "T0var21", "T0var22","m1", "m2")

oscillatingm <- ctModel(n.latent = 2, n.manifest = 1, Tpoints = 11, 
  MANIFESTVAR = matrix(c(0), nrow = 1, ncol = 1), 
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  T0MEANS = matrix(c('m1', 'm2'), nrow = 2, ncol = 1), 
  T0VAR = matrix(c("T0var11", "T0var21", 0, "T0var22"), nrow = 2, ncol = 2),
  DRIFT = matrix(c(0, "crosseffect", 1, "autoeffect"), nrow = 2, ncol = 2), 
  CINT = matrix(0, ncol = 1, nrow = 2),
  DIFFUSION = matrix(c(0, 0, 0, "diffusion"), nrow = 2, ncol = 2))

oscillatingf <- ctFit(Oscillating, oscillatingm, carefulFit = FALSE,retryattempts=1)
mxs <- summary(oscillatingf)

sm <- ctStanModel(oscillatingm)
sm$pars$indvarying <- FALSE

ld <- ctWideToLong(Oscillating,Tpoints = 11,n.manifest = 1,n.TDpred = 0,n.TIpred = 0)
ld <- ctDeintervalise(ld)

sf <- ctStanFit(datalong = ld, ctstanmodel = sm,optimize=TRUE,verbose=1,optimcontrol=list(isloops=0,deoptim=FALSE),
  nlcontrol=list(nldynamics=FALSE),nopriors = T,chains = 2,cores=2)
s=summary(sf)

expect_equal(-3461.936,oscillatingf$mxobj$output$Minus2LogLikelihood,tolerance=.001)
expect_equivalent(s$popmeans[order(rownames(s$popmeans)),'mean'],
  mxs$ctparameters[order(rownames(mxs$ctparameters)),'Value'],
  tolerance=.05)


})
}
