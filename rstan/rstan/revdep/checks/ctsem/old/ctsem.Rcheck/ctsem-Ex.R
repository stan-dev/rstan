pkgname <- "ctsem"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ctsem')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Kalman")
### * Kalman

flush(stderr()); flush(stdout())

### Name: Kalman
### Title: Kalman
### Aliases: Kalman

### ** Examples

### ctstantestfit is a dummy ctStanFit object with 2 manifest indicators,
###  4 latents, and 1 time dependent predictor.

### get parameter matrices
kpars <- ctStanContinuousPars(ctstantestfit)

#construct dummy data
datalong <- cbind(0:9, 1, matrix(rnorm(20,2,1),ncol=2))
datalong[c(1:3,9:10),3:4]<-NA #missing data to pre/fore cast
colnames(datalong) <- c('time', 'id', paste0('Y',1:2))
print(datalong)

#obtain Kalman filtered estimates
kout <- Kalman(kpars=kpars, datalong=datalong,
  manifestNames=paste0('Y',1:nrow(kpars$MANIFESTMEANS)),
  latentNames=paste0('eta',1:nrow(kpars$DRIFT)))

#print and plot smoothed estimates (conditional on all states) of indicators.
print(kout$ysmooth)
matplot(kout$time,kout$ysmooth,type='l')
matplot(kout$time,datalong[,3:4],type='p',add=TRUE,pch=1)



cleanEx()
nameEx("ctCI")
### * ctCI

flush(stderr()); flush(stdout())

### Name: ctCI
### Title: ctCI Computes confidence intervals on specified parameters /
###   matrices for already fitted ctsem fit object.
### Aliases: ctCI

### ** Examples

## Examples set to 'dontrun' because they take longer than 5s.
## Not run: 
##D data("ctExample3")
##D model <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 100, 
##D  LAMBDA = matrix(c(1, "lambda2", "lambda3"), nrow = 3, ncol = 1), 
##D  MANIFESTMEANS = matrix(c(0, "manifestmean2", "manifestmean3"), nrow = 3, 
##D    ncol = 1))
##D fit <- ctFit(dat = ctExample3, ctmodelobj = model, objective = "Kalman",
##D  stationary = c("T0VAR"))
##D 
##D fit <- ctCI(fit, confidenceintervals = 'DRIFT')
##D 
##D summary(fit)$omxsummary$CI
## End(Not run)



cleanEx()
nameEx("ctCheckFit")
### * ctCheckFit

flush(stderr()); flush(stdout())

### Name: ctCheckFit
### Title: Check absolute fit of ctFit or ctStanFit object.
### Aliases: ctCheckFit

### ** Examples

## Not run: 
##D data(ctExample1)
##D traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
##D   manifestNames=c('LeisureTime', 'Happiness'), 
##D   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
##D traitfit <- ctFit(dat=ctExample1, ctmodelobj=traitmodel)
##D 
##D check <- ctCheckFit(traitfit,niter=5)
##D plot(check)
## End(Not run)



cleanEx()
nameEx("ctCollapse")
### * ctCollapse

flush(stderr()); flush(stdout())

### Name: ctCollapse
### Title: ctCollapse Easily collapse an array margin using a specified
###   function.
### Aliases: ctCollapse

### ** Examples

testarray <- array(rnorm(900,2,1),dim=c(100,3,3))
ctCollapse(testarray,1,mean)



cleanEx()
nameEx("ctDensity")
### * ctDensity

flush(stderr()); flush(stdout())

### Name: ctDensity
### Title: ctDensity
### Aliases: ctDensity

### ** Examples

y <- ctDensity(exp(rnorm(80)))
plot(y$density,xlim=y$xlim,ylim=y$ylim)

#### Compare to base defaults:
par(mfrow=c(1,2))
y=exp(rnorm(10000))
ctdens<-ctDensity(y)
plot(ctdens$density, ylim=ctdens$ylim,xlim=ctdens$xlim)
plot(density(y))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ctDiscretePars")
### * ctDiscretePars

flush(stderr()); flush(stdout())

### Name: ctDiscretePars
### Title: ctDiscretePars
### Aliases: ctDiscretePars

### ** Examples

pars <- ctStanContinuousPars(ctstantestfit)
ctDiscretePars(pars,times=c(.5,1))




cleanEx()
nameEx("ctDiscretiseData")
### * ctDiscretiseData

flush(stderr()); flush(stdout())

### Name: ctDiscretiseData
### Title: Discretise long format continuous time (ctsem) data to specific
###   timestep.
### Aliases: ctDiscretiseData

### ** Examples

long <- ctWideToLong(datawide=ctExample2,Tpoints=8,n.manifest=2,n.TDpred=1,
 manifestNames=c('LeisureTime','Happiness'),
 TDpredNames=c('MoneyInt'))

long <- ctDeintervalise(long)

long <- ctDiscretiseData(dlong=long, timestep = 1.1,TDpredNames=c('MoneyInt'))



cleanEx()
nameEx("ctExample2")
### * ctExample2

flush(stderr()); flush(stdout())

### Name: ctExample2
### Title: ctExample2
### Aliases: ctExample2

### ** Examples

## Not run: 
##D #two process, one time dependent predictor example
##D Tpoints=20
##D manifestNames<-c('LeisureTime','Happiness')
##D TDpredNames<-'MoneyInt'
##D testm<-ctModel(Tpoints=Tpoints,n.latent=3,n.TDpred=1,n.TIpred=0,n.manifest=2,    
##D   LAMBDA=cbind(diag(1,2),0),
##D   MANIFESTVAR=diag(.1,2),
##D   DRIFT=matrix(c(-.3,.12,0,  -.02,-.3,0, 1,-.3,-.0001  ),nrow=3,ncol=3),
##D   TRAITVAR=t(chol(matrix(c(.2,-.1,0,  -.1,.21,0,  0,0,0.00001),ncol=3,nrow=3))),
##D   DIFFUSION=t(chol(diag(c(1.2,.6,0.0001),3))),
##D   CINT=matrix(c(1,.3,0),nrow=3),
##D   T0MEANS=matrix(0,ncol=1,nrow=3),
##D   T0VAR=diag(c(1,1,0),3),
##D   TDPREDEFFECT=matrix(c(.6,.4,1),nrow=3),
##D   TDPREDVAR=diag(c(rep(0,Tpoints)),Tpoints),
##D   TDPREDMEANS=matrix(c(0,0,0,0,0,1,rep(0,Tpoints-6)),ncol=1,nrow=(Tpoints)))
##D testd<-ctGenerate(testm,n.subjects=10,burnin=10) #generate data
##D 
##D ctIndplot(testd,Tpoints=Tpoints,n.manifest=2,n.subjects=10,colourby="variable")
##D 
##D timestokeep=c(0,1,4,5,7,8,16,19)
##D deltaT<-timestokeep[-1] - timestokeep[-8]
##D testd<-testd[,c(paste0('Y',1:2,'_T',rep(timestokeep,each=2)),paste0('TD1_T',timestokeep))]
##D testd<-cbind(testd,matrix(deltaT,nrow=nrow(testd),ncol=length(deltaT),byrow=TRUE))
##D 
##D colnames(testd)<-ctWideNames(n.manifest=2,Tpoints=8,n.TDpred=1,
##D manifestNames=manifestNames,TDpredNames=TDpredNames)
##D ctExample2<-testd
##D save(ctExample2,file=".\\data\\ctExample2.rda") 
## End(Not run)



cleanEx()
nameEx("ctExample4")
### * ctExample4

flush(stderr()); flush(stdout())

### Name: ctExample4
### Title: ctExample4
### Aliases: ctExample4

### ** Examples

## Not run: 
##D Tpoints=20
##D subjects=20
##D full<-c()
##D for(i in 1:20){
##D   LAMBDA<-matrix(c(1,.7,ifelse(i >(subjects/2),.2,1.4)))
##D   print(LAMBDA)
##D   testm<-ctModel(Tpoints=Tpoints,n.latent=1,n.TDpred=0,n.TIpred=0,n.manifest=3,
##D     MANIFESTVAR=diag(.2,3),
##D     # TRAITVAR=diag(.2,1),
##D     LAMBDA=LAMBDA,
##D     DRIFT=matrix(c(-.1),nrow=1,ncol=1),
##D     DIFFUSION=diag(c(.12),1),
##D     MANIFESTMEANS=matrix(c(0,.42,1.3),ncol=1),
##D     CINT=matrix(c(.2),ncol=1))
##D   
##D   testd<-ctGenerate(testm,n.subjects=1,burnin=300)
##D   full<-rbind(full,testd)
##D }
##D 
##D ctExample4<-full
##D save(ctExample4,file=".\\data\\ctExample4.rda") #save wide format example
## End(Not run)



cleanEx()
nameEx("ctFit")
### * ctFit

flush(stderr()); flush(stdout())

### Name: ctFit
### Title: Fit a ctsem object
### Aliases: ctFit

### ** Examples

## Examples set to 'dontrun' because they take longer than 5s.
## Not run: 
##D mfrowOld<-par()$mfrow
##D par(mfrow=c(2, 3))
##D 
##D ### example from Driver, Oud, Voelkle (2017), 
##D ### simulated happiness and leisure time with unobserved heterogeneity.
##D data(ctExample1)
##D traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
##D   manifestNames=c('LeisureTime', 'Happiness'), 
##D   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
##D traitfit <- ctFit(dat=ctExample1, ctmodelobj=traitmodel)
##D summary(traitfit)
##D plot(traitfit, wait=FALSE)
##D 
##D 
##D ###Example from Voelkle, Oud, Davidov, and Schmidt (2012) - anomia and authoritarianism.  
##D data(AnomAuth) 
##D AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
##D Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = NULL) 
##D AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)
##D summary(AnomAuthfit)
##D 
##D 
##D ### Single subject time series - using Kalman filter (OpenMx statespace expectation)
##D data('ctExample3')
##D model <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 100, 
##D   LAMBDA = matrix(c(1, 'lambda2', 'lambda3'), nrow = 3, ncol = 1), 
##D   CINT= matrix('cint'),
##D   MANIFESTMEANS = matrix(c(0, 'manifestmean2', 'manifestmean3'), nrow = 3, 
##D     ncol = 1))
##D fit <- ctFit(dat = ctExample3, ctmodelobj = model, objective = 'Kalman', 
##D   stationary = c('T0VAR'))
##D 
##D 
##D ###Oscillating model from Voelkle & Oud (2013). 
##D data("Oscillating")
##D 
##D inits <- c(-39, -.3, 1.01, 10.01, .1, 10.01, 0.05, .9, 0)
##D names(inits) <- c("crosseffect","autoeffect", "diffusion",
##D   "T0var11", "T0var21", "T0var22","m1", "m2", 'manifestmean')
##D 
##D oscillatingm <- ctModel(n.latent = 2, n.manifest = 1, Tpoints = 11,
##D   MANIFESTVAR = matrix(c(0), nrow = 1, ncol = 1),
##D   LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
##D   T0MEANS = matrix(c('m1', 'm2'), nrow = 2, ncol = 1),
##D   T0VAR = matrix(c("T0var11", "T0var21", 0, "T0var22"), nrow = 2, ncol = 2),
##D   DRIFT = matrix(c(0, "crosseffect", 1, "autoeffect"), nrow = 2, ncol = 2),
##D   CINT = matrix(0, ncol = 1, nrow = 2),
##D   MANIFESTMEANS = matrix('manifestmean', nrow = 1, ncol = 1),
##D   DIFFUSION = matrix(c(0, 0, 0, "diffusion"), nrow = 2, ncol = 2),
##D   startValues=inits)
##D 
##D oscillatingf <- ctFit(Oscillating, oscillatingm, carefulFit = FALSE)
## End(Not run)



cleanEx()
nameEx("ctGenerate")
### * ctGenerate

flush(stderr()); flush(stdout())

### Name: ctGenerate
### Title: ctGenerate
### Aliases: ctGenerate

### ** Examples

#generate data for 2 process model, each process measured by noisy indicator, 
#stable individual differences in process levels.

generatingModel<-ctModel(Tpoints=8,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
 MANIFESTVAR=diag(.1,2),
 LAMBDA=diag(1,2),
 DRIFT=matrix(c(-.2,-.05,-.1,-.1),nrow=2),
 TRAITVAR=matrix(c(.5,.2,0,.8),nrow=2),
 DIFFUSION=matrix(c(1,.2,0,4),2),
 CINT=matrix(c(1,0),nrow=2),
 T0MEANS=matrix(0,ncol=1,nrow=2),
 T0VAR=diag(1,2))

data<-ctGenerate(generatingModel,n.subjects=15,burnin=10)



cleanEx()
nameEx("ctGenerateFromFit")
### * ctGenerateFromFit

flush(stderr()); flush(stdout())

### Name: ctGenerateFromFit
### Title: Generates data according to the model estimated in a ctsemFit
###   object.#'
### Aliases: ctGenerateFromFit

### ** Examples


data(AnomAuth) 
AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
  Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2)) 
AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)

dwide <- ctGenerateFromFit(AnomAuthfit,timestep=1,n.subjects=5)

par(mfrow=c(1,2))
ctIndplot(datawide = dwide,n.subjects = 5,n.manifest = 2,vars=1,Tpoints = 4)
ctIndplot(datawide = AnomAuth+rnorm(length(AnomAuth)),vars=1,n.subjects = 5,
n.manifest = 2,Tpoints = 4)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ctIndplot")
### * ctIndplot

flush(stderr()); flush(stdout())

### Name: ctIndplot
### Title: ctIndplot
### Aliases: ctIndplot

### ** Examples


data(ctExample1)
ctIndplot(ctExample1,n.subjects=1, n.manifest=2,Tpoints=6, colourby='variable')




cleanEx()
nameEx("ctIntervalise")
### * ctIntervalise

flush(stderr()); flush(stdout())

### Name: ctIntervalise
### Title: Converts absolute times to intervals for wide format ctsem panel
###   data
### Aliases: ctIntervalise

### ** Examples

#First load the long format data with absolute times
data('longexample')

#Then convert to wide format
wideexample <- ctLongToWide(datalong = longexample, id = "id", 
time = "time", manifestNames = c("Y1", "Y2", "Y3"), 
TDpredNames = "TD1", TIpredNames = c("TI1", "TI2"))

#Then convert the absolute times to intervals, using the Tpoints reported from the prior step.
wide <- ctIntervalise(datawide = wideexample, Tpoints = 4, n.manifest = 3, 
n.TDpred = 1, n.TIpred = 2, manifestNames = c("Y1", "Y2", "Y3"), 
TDpredNames = "TD1", TIpredNames = c("TI1", "TI2") )
 
print(wide)



cleanEx()
nameEx("ctKalman")
### * ctKalman

flush(stderr()); flush(stdout())

### Name: ctKalman
### Title: ctKalman
### Aliases: ctKalman

### ** Examples

## Not run: 
##D #Basic
##D ctKalman(ctstantestfit, timerange=c(0,60), timestep=.5, plot=TRUE)
##D 
##D #Multiple subjects, y and yprior, showing plot arguments
##D ctKalman(ctstantestfit, timerange=c(0,60), timestep=.1, plot=TRUE,
##D   subjects=2:3, 
##D   kalmanvec=c('y','yprior'),
##D   errorvec=c(NA,'ypriorcov'), #'auto' would also have achieved this
##D   ltyvec="auto",
##D   colvec='auto', 
##D   lwdvec='auto', 
##D   subsetindices=2, #Only plotting 2nd dimension of y and yprior
##D   pchvec='auto', typevec='auto',grid=TRUE,legend=TRUE,
##D   plotcontrol=list(xlim=c(0,55),main='Observations and priors'),
##D   polygoncontrol=list(steps=5))
##D   
## End(Not run)



cleanEx()
nameEx("ctKalmanPlot")
### * ctKalmanPlot

flush(stderr()); flush(stdout())

### Name: ctKalmanPlot
### Title: ctKalmanPlot
### Aliases: ctKalmanPlot

### ** Examples

### Get output from ctKalman
x<-ctKalman(ctstantestfit,subjects=2)

### Plot with ctKalmanPlot
ctKalmanPlot(x, subjects=2)

###Single step procedure:
ctKalman(ctstantestfit,subjects=2,plot=TRUE)



cleanEx()
nameEx("ctLongToWide")
### * ctLongToWide

flush(stderr()); flush(stdout())

### Name: ctLongToWide
### Title: ctLongToWide Restructures time series / panel data from long
###   format to wide format for ctsem analysis
### Aliases: ctLongToWide

### ** Examples

#First load the long format data with absolute times
data('longexample')

#Then convert to wide format
wideexample <- ctLongToWide(datalong = longexample, id = "id", 
time = "time", manifestNames = c("Y1", "Y2", "Y3"), 
TDpredNames = "TD1", TIpredNames = c("TI1", "TI2"))

#Then convert the absolute times to intervals, using the Tpoints reported from the prior step.
wide <- ctIntervalise(datawide = wideexample, Tpoints = 4, n.manifest = 3, 
n.TDpred = 1, n.TIpred = 2, manifestNames = c("Y1", "Y2", "Y3"), 
TDpredNames = "TD1", TIpredNames = c("TI1", "TI2") )




cleanEx()
nameEx("ctModel")
### * ctModel

flush(stderr()); flush(stdout())

### Name: ctModel
### Title: Define a ctsem model
### Aliases: ctModel

### ** Examples

 ### Frequentist example:
 ### impulse and level change time dependent predictor 
 ### example from Driver, Oud, Voelkle (2015)
 data('ctExample2')
 tdpredmodel <- ctModel(n.manifest = 2, n.latent = 3, n.TDpred = 1, 
 Tpoints = 8, manifestNames = c('LeisureTime', 'Happiness'), 
 TDpredNames = 'MoneyInt', 
 latentNames = c('LeisureTime', 'Happiness', 'MoneyIntLatent'),
 LAMBDA = matrix(c(1,0, 0,1, 0,0), ncol = 3), TRAITVAR = "auto")

 tdpredmodel$TRAITVAR[3, ] <- 0
 tdpredmodel$TRAITVAR[, 3] <- 0
 tdpredmodel$DIFFUSION[, 3] <- 0
 tdpredmodel$DIFFUSION[3, ] <- 0
 tdpredmodel$T0VAR[3, ] <- 0
 tdpredmodel$T0VAR[, 3] <- 0
 tdpredmodel$CINT[3] <- 0
 tdpredmodel$T0MEANS[3] <- 0
 tdpredmodel$TDPREDEFFECT[3, ] <- 1
 tdpredmodel$DRIFT[3, ] <- 0
 
 
###Bayesian example:
model<-ctModel(type='stanct',
n.latent=2, latentNames=c('eta1','eta2'),
n.manifest=2, manifestNames=c('Y1','Y2'),
n.TDpred=1, TDpredNames='TD1', 
n.TIpred=3, TIpredNames=c('TI1','TI2','TI3'),
LAMBDA=diag(2))





cleanEx()
nameEx("ctModelFromFit")
### * ctModelFromFit

flush(stderr()); flush(stdout())

### Name: ctModelFromFit
### Title: Extract a ctsem model structure with parameter values from a
###   ctsem fit object.
### Aliases: ctModelFromFit

### ** Examples

data(AnomAuth) 
AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
  Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2)) 
AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)

fitmodel <- ctModelFromFit(AnomAuthfit)



cleanEx()
nameEx("ctMultigroupFit")
### * ctMultigroupFit

flush(stderr()); flush(stdout())

### Name: ctMultigroupFit
### Title: Fits a multiple group continuous time model.
### Aliases: ctMultigroupFit

### ** Examples

## Not run: 
##D 
##D #Two group model, all parameters except LAMBDA[3,1] constrained across groups.
##D data(ctExample4)
##D basemodel<-ctModel(n.latent=1, n.manifest=3, Tpoints=20,
##D                    LAMBDA=matrix(c(1, 'lambda2', 'lambda3'), nrow=3, ncol=1),
##D                    MANIFESTMEANS=matrix(c(0, 'manifestmean2', 'manifestmean3'), 
##D                    nrow=3, ncol=1), TRAITVAR = 'auto')
##D 
##D freemodel<-basemodel
##D freemodel$LAMBDA[3,1]<-'groupfree'
##D groups<-paste0('g',rep(1:2, each=10),'_')
##D 
##D multif<-ctMultigroupFit(dat=ctExample4, groupings=groups,
##D                        ctmodelobj=basemodel, freemodel=freemodel)
##D summary(multif,group=1)
##D 
##D 
##D 
##D #fixed model approach
##D fixedmodel<-basemodel
##D fixedmodel$LAMBDA[2,1]<-'groupfixed'
##D groups<-paste0('g',rep(1:2, each=10),'_')
##D 
##D multif<-ctMultigroupFit(dat=ctExample4, groupings=groups,
##D                        ctmodelobj=basemodel, fixedmodel=fixedmodel)
##D summary(multif) 
## End(Not run)





cleanEx()
nameEx("ctPlot")
### * ctPlot

flush(stderr()); flush(stdout())

### Name: ctPlot
### Title: ctPlot
### Aliases: ctPlot

### ** Examples

## Examples set to 'dontrun' because they take longer than 5s.

### example from Driver, Oud, Voelkle (2016), 
### simulated happiness and leisure time with unobserved heterogeneity.
## Not run: 
##D data(ctExample1)
##D traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
##D   manifestNames=c('LeisureTime', 'Happiness'), 
##D   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
##D traitfit <- ctFit(dat=ctExample1, ctmodelobj=traitmodel)
##D ctPlot(traitfit, plotType='CR', xlim=c(0,5),ylim=c(-1,1))
## End(Not run)



cleanEx()
nameEx("ctPlotArray")
### * ctPlotArray

flush(stderr()); flush(stdout())

### Name: ctPlotArray
### Title: Plots three dimensional y values for quantile plots
### Aliases: ctPlotArray

### ** Examples

## Not run: 
##D input<-ctStanTIpredeffects(ctstantestfit, plot=FALSE, whichpars='CINT', 
##D  nsamples=10,nsubjects=10)
##D     
##D ctPlotArray(input=input)
## End(Not run)



cleanEx()
nameEx("ctPoly")
### * ctPoly

flush(stderr()); flush(stdout())

### Name: ctPoly
### Title: Plots uncertainty bands with shading
### Aliases: ctPoly

### ** Examples

plot(0:100,sqrt(0:100),type='l')
ctPoly(x=0:100, y=sqrt(0:100), 
yhigh=sqrt(0:100) - runif(101), 
ylow=sqrt(0:100) + runif(101),
col=adjustcolor('red',alpha.f=.1))



cleanEx()
nameEx("ctPostPredict")
### * ctPostPredict

flush(stderr()); flush(stdout())

### Name: ctPostPredict
### Title: Posterior predictive type check for ctsemFit.
### Aliases: ctPostPredict

### ** Examples

data("AnomAuth")
AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
  Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = 'auto') 
AnomAuthFit <- ctFit(AnomAuth, AnomAuthmodel)
ctPostPredict(AnomAuthFit,timestep=.5,n.subjects=100)



cleanEx()
nameEx("ctStanContinuousPars")
### * ctStanContinuousPars

flush(stderr()); flush(stdout())

### Name: ctStanContinuousPars
### Title: ctStanContinuousPars
### Aliases: ctStanContinuousPars

### ** Examples

## Not run: 
##D #posterior median over all subjects (also reflects mean of unconstrained pars)
##D ctStanContinuousPars(ctstantestfit)
##D 
##D #posterior 97.5% quantiles for subject 2
##D ctStanContinuousPars(ctstantestfit, subjects=2, calcfunc=quantile, 
##D calcfuncargs=list(probs=0.975))
## End(Not run)



cleanEx()
nameEx("ctStanDiscretePars")
### * ctStanDiscretePars

flush(stderr()); flush(stdout())

### Name: ctStanDiscretePars
### Title: ctStanDiscretePars
### Aliases: ctStanDiscretePars

### ** Examples

ctStanDiscretePars(ctstantestfit,times=seq(.5,4,.1), 
 plot=TRUE,indices='all')



cleanEx()
nameEx("ctStanDiscreteParsPlot")
### * ctStanDiscreteParsPlot

flush(stderr()); flush(stdout())

### Name: ctStanDiscreteParsPlot
### Title: ctStanDiscreteParsPlot
### Aliases: ctStanDiscreteParsPlot

### ** Examples

x <- ctStanDiscretePars(ctstantestfit)

ctStanDiscreteParsPlot(x, 'CR')



cleanEx()
nameEx("ctStanFit")
### * ctStanFit

flush(stderr()); flush(stdout())

### Name: ctStanFit
### Title: ctStanFit
### Aliases: ctStanFit

### ** Examples

## Not run: 
##D #test data with 2 manifest indicators measuring 1 latent process each, 
##D # 1 time dependent predictor, 3 time independent predictors
##D head(ctstantestdat) 
##D 
##D #generate a ctStanModel
##D model<-ctModel(type='stanct',
##D n.latent=2, latentNames=c('eta1','eta2'),
##D n.manifest=2, manifestNames=c('Y1','Y2'),
##D n.TDpred=1, TDpredNames='TD1', 
##D n.TIpred=3, TIpredNames=c('TI1','TI2','TI3'),
##D LAMBDA=diag(2))
##D 
##D #set all parameters except manifest means to be fixed across subjects
##D model$pars$indvarying[-c(19,20)] <- FALSE
##D 
##D #fit model to data (takes a few minutes - but insufficient 
##D # iterations and max_treedepth for inference!)
##D fit<-ctStanFit(ctstantestdat, model, iter=200, chains=2, 
##D control=list(max_treedepth=6))
##D 
##D #output functions
##D summary(fit) 
##D 
##D plot(fit)
##D 
##D 
##D 
##D ###### EXTENDED EXAMPLES #######
##D 
##D 
##D 
##D library(ctsem)
##D set.seed(3)
##D 
##D #recommended to adjust these to appropriate number of cores on machine / 
##D # chains desired. (min 3 chains recommended, but not necessary here)
##D setcores <- 3
##D setchains <- 3
##D 
##D #### Data generation (this section needs to be run, but not necessary to understand!)
##D Tpoints <- 20
##D nmanifest <- 4
##D nlatent <- 2
##D nsubjects<-20
##D 
##D #random effects
##D age <- rnorm(nsubjects) #standardised
##D cint1<-rnorm(nsubjects,2,.3)+age*.5
##D cint2 <- cint1*.5+rnorm(nsubjects,1,.2)+age*.5
##D tdpredeffect <- rnorm(nsubjects,5,.3)+age*.5
##D 
##D for(i in 1:nsubjects){
##D   #generating model
##D   gm<-ctModel(Tpoints=Tpoints,n.manifest = nmanifest,n.latent = nlatent,n.TDpred = 1,
##D     LAMBDA = matrix(c(1,0,0,0, 0,1,.8,1.3),nrow=nmanifest,ncol=nlatent),
##D     DRIFT=matrix(c(-.3, .1, 0, -.5),nlatent,nlatent),
##D     TDPREDEFFECT=matrix(c(tdpredeffect[i],0),nrow=nlatent),
##D     TDPREDMEANS=matrix(c(rep(0,Tpoints-10),1,rep(0,9)),ncol=1),
##D     DIFFUSION = matrix(c(.5, 0, 0, .5),2,2),
##D     CINT = matrix(c(cint1[i],cint2[i]),ncol=1),
##D     T0VAR=diag(2,nlatent,nlatent),
##D     MANIFESTVAR = diag(.5, nmanifest))
##D   
##D   #generate data
##D   newdat <- ctGenerate(ctmodelobj = gm,n.subjects = 1,burnin = 2,
##D     dtmat<-rbind(c(rep(.5,8),3,rep(.5,Tpoints-9))),
##D       wide = FALSE)
##D   newdat[,'id'] <- i #set id for each subject
##D   newdat <- cbind(newdat,age[i]) #include time independent predictor
##D   if(i ==1) {
##D     dat <- newdat[1:(Tpoints-10),] #pre intervention data
##D     dat2 <- newdat #including post intervention data
##D   }
##D   if(i > 1) {
##D     dat <- rbind(dat, newdat[1:(Tpoints-10),])
##D     dat2 <- rbind(dat2,newdat)
##D   }
##D }
##D colnames(dat)[ncol(dat)] <- 'age'
##D colnames(dat2)[ncol(dat)] <- 'age'
##D 
##D 
##D #plot generated data for sanity
##D plot(age)
##D matplot(dat[,gm$manifestNames],type='l',pch=1)
##D plotvar <- 'Y1'
##D plot(dat[dat[,'id']==1,'time'],dat[dat[,'id']==1,plotvar],type='l',
##D   ylim=range(dat[,plotvar],na.rm=TRUE))
##D for(i in 2:nsubjects){
##D   points(dat[dat[,'id']==i,'time'],dat[dat[,'id']==i,plotvar],type='l',col=i)
##D }
##D 
##D 
##D 
##D 
##D 
##D #### Model fitting (from here it is good to understand!)
##D 
##D #Specify univariate linear growth curve
##D #page 5 of https://cran.r-project.org/web/packages/ctsem/vignettes/hierarchical.pdf 
##D # documents these arguments (or use ?ctModel )
##D 
##D m1 <- ctModel(n.manifest = 1,n.latent = 1,n.TIpred = 1, type = 'stanct',
##D   manifestNames = c('Y1'), latentNames=c('L1'),TIpredNames = 'age',
##D   DRIFT=matrix(-1e-5,nrow=1,ncol=1),
##D   DIFFUSION=matrix(0,nrow=1,ncol=1),
##D   CINT=matrix(c('cint1'),ncol=1),
##D   T0MEANS=matrix(c('t0m1'),ncol=1),
##D   T0VAR=matrix(0,nrow=1,ncol=1),
##D   LAMBDA = diag(1),
##D   MANIFESTMEANS=matrix(0,ncol=1),
##D   MANIFESTVAR=matrix(c('merror1'),nrow=1,ncol=1))
##D 
##D #modify between subject aspects -- alternatively, run: edit(m1$pars)
##D m1$pars$indvarying[-which(m1$pars$matrix %in% c('T0MEANS','CINT'))] <- FALSE
##D m1$pars$age_effect[-which(m1$pars$matrix %in% c('T0MEANS','CINT'))] <- FALSE
##D 
##D 
##D plot(m1) #plot prior distributions
##D 
##D #fit
##D f1 <- ctStanFit(datalong = dat, ctstanmodel = m1, 
##D   cores = setcores,chains = setchains,plot=TRUE,
##D   control=list(max_treedepth=7),iter=150)
##D 
##D summary(f1)
##D 
##D #plots of individual subject models v data
##D ctKalman(f1,timestep=.01,plot=TRUE,subjects=1:4,kalmanvec=c('y','etasmooth'))
##D ctKalman(f1,timestep=.01,plot=TRUE,subjects=1,kalmanvec=c('y','ysmooth'))
##D 
##D ctStanPlotPost(f1) #compare prior to posterior distributions 
##D ctStanPlotPost(f1, priorwidth = FALSE) #rescale to width of posterior 
##D 
##D ctStanPostPredict(f1) #compare randomly generated data from posterior to observed data
##D 
##D cf<-ctCheckFit(f1) #compare covariance of randomly generated data to observed cov
##D plot(cf)
##D 
##D #accessing the stan object directly 
##D library(rstan)
##D postsamples <- extract(f1$stanfit,pars='Ygen') #extract data generated from posterior
##D plot( f1$data$time, 
##D   postsamples$Ygen[1,1, ,1]) #iteration 1 (already shuffled), chain 1, all occasions, var 1.
##D points(f1$data$time, f1$data$Y[,1],col='red') #1st manifest variable
##D 
##D 
##D 
##D 
##D 
##D #Specify model including dynamics
##D m2 <- ctModel(n.manifest = 1,n.latent = 1,n.TIpred = 1, type = 'stanct',
##D   manifestNames = c('Y1'), latentNames=c('L1'),TIpredNames = 'age',
##D   DRIFT=matrix('drift11',nrow=1,ncol=1),
##D   DIFFUSION=matrix('diffusion11',nrow=1,ncol=1),
##D   CINT=matrix(c('cint1'),ncol=1),
##D   T0MEANS=matrix(c('t0m1'),ncol=1),
##D   T0VAR=matrix('t0var11',nrow=1,ncol=1),
##D   LAMBDA = diag(1),
##D   MANIFESTMEANS=matrix(0,ncol=1),
##D   MANIFESTVAR=matrix(c('merror1'),nrow=1,ncol=1))
##D 
##D m2$pars$indvarying[-which(m2$pars$matrix %in% c('T0MEANS','CINT'))] <- FALSE
##D m2$pars$age_effect[-which(m2$pars$matrix %in% c('T0MEANS','CINT'))] <- FALSE
##D 
##D 
##D f2 <- ctStanFit(datalong = dat, ctstanmodel = m2, cores = setcores,
##D   chains = setchains,plot=TRUE,
##D   control=list(max_treedepth=7),iter=150)
##D 
##D summary(f2,parmatrices=TRUE,timeinterval=1)
##D 
##D ctKalman(f2,timestep=.01,plot=TRUE,subjects=1,kalmanvec=c('y','etaprior'))
##D ctKalman(f2,timestep=.01,plot=TRUE,subjects=1:4,kalmanvec=c('y','etasmooth'))
##D ctKalman(f2,timestep=.01,plot=TRUE,subjects=1:2,kalmanvec=c('y','ysmooth'))
##D 
##D ctStanPlotPost(f2)
##D 
##D ctStanPostPredict(f2)
##D 
##D 
##D 
##D 
##D 
##D #Include intervention
##D m3 <- ctModel(n.manifest = 1,n.latent = 1,n.TIpred = 1, type = 'stanct',
##D   manifestNames = c('Y1'), latentNames=c('L1'),TIpredNames = 'age',
##D   n.TDpred=1,TDpredNames = 'TD1', #this line includes the intervention
##D   TDPREDEFFECT=matrix(c('tdpredeffect'),nrow=1,ncol=1), #intervention effect
##D   DRIFT=matrix('drift11',nrow=1,ncol=1),
##D   DIFFUSION=matrix('diffusion11',nrow=1,ncol=1),
##D   CINT=matrix(c('cint1'),ncol=1),
##D   T0MEANS=matrix(c('t0m1'),ncol=1),
##D   T0VAR=matrix('t0var11',nrow=1,ncol=1),
##D   LAMBDA = diag(1),
##D   MANIFESTMEANS=matrix(0,ncol=1),
##D   MANIFESTVAR=matrix(c('merror1'),nrow=1,ncol=1))
##D 
##D m3$pars$indvarying[-which(m3$pars$matrix %in% 
##D   c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE
##D   
##D m3$pars$age_effect[-which(m3$pars$matrix %in% 
##D   c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE
##D 
##D f3 <- ctStanFit(datalong = dat2, ctstanmodel = m3, cores = setcores,
##D   chains = setchains,plot=TRUE,
##D   control=list(max_treedepth=7),iter=150)
##D 
##D summary(f3,parmatrices=TRUE)
##D 
##D ctKalman(f3,timestep=.01,plot=TRUE,subjects=1,kalmanvec=c('y','etaprior'))
##D ctKalman(f3,timestep=.01,plot=TRUE,subjects=1:4,kalmanvec=c('y','etasmooth'))
##D ctKalman(f3,timestep=.01,plot=TRUE,subjects=1:2,kalmanvec=c('y','ysmooth'))
##D 
##D ctStanPlotPost(f3)
##D 
##D ctStanPostPredict(f3, datarows=0:100)
##D 
##D 
##D 
##D 
##D 
##D 
##D #include 2nd latent process
##D 
##D #use either full explicit specification
##D m4 <- ctModel(n.manifest = 2,n.latent = 2,n.TIpred = 1, type = 'stanct', #no of vars updated
##D   manifestNames = c('Y1','Y2'), latentNames=c('L1','L2'),TIpredNames = 'age', 
##D   n.TDpred=1,TDpredNames = 'TD1', 
##D   TDPREDEFFECT=matrix(c('tdpredeffect1','tdpredeffect2'),nrow=2,ncol=1), 
##D   DRIFT=matrix(c('drift11','drift21','drift12','drift22'),nrow=2,ncol=2),
##D   DIFFUSION=matrix(c('diffusion11','diffusion21',0,'diffusion22'),nrow=2,ncol=2),
##D   CINT=matrix(c('cint1','cint2'),nrow=2,ncol=1),
##D   T0MEANS=matrix(c('t0m1','t0m2'),nrow=2,ncol=1),
##D   T0VAR=matrix(c('t0var11','t0var21',0,'t0var22'),nrow=2,ncol=2),
##D   LAMBDA = matrix(c(1,0,0,1),nrow=2,ncol=2),
##D   MANIFESTMEANS=matrix(c(0,0),nrow=2,ncol=1),
##D   MANIFESTVAR=matrix(c('merror1',0,0,'merror2'),nrow=2,ncol=2))
##D 
##D #restrict between subjects variation / covariate effects
##D m4$pars$indvarying[-which(m4$pars$matrix %in% c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE
##D m4$pars$age_effect[-which(m4$pars$matrix %in% c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE
##D 
##D #or rely on defaults (MANIFESTMEANS now free instead of CINT -- 
##D #  no substantive difference for one indicator factors)
##D m4 <- ctModel(n.manifest = 2,n.latent = 2,n.TIpred = 1, type = 'stanct', 
##D   manifestNames = c('Y1','Y2'), latentNames=c('L1','L2'),TIpredNames = 'age',
##D   n.TDpred=1,TDpredNames = 'TD1',
##D   LAMBDA = matrix(c(1,0,0,1),nrow=2,ncol=2))
##D 
##D #restrict between subjects variation / covariate effects
##D m4$pars$indvarying[-which(m4$pars$matrix %in% c('T0MEANS','MANIFESTMEANS','TDPREDEFFECT'))] <- FALSE
##D m4$pars$age_effect[-which(m4$pars$matrix %in% c('T0MEANS','MANIFESTMEANS','TDPREDEFFECT'))] <- FALSE
##D 
##D f4 <- ctStanFit(datalong = dat2, ctstanmodel = m4, cores = setcores,chains = setchains,plot=TRUE,
##D   optimize=TRUE,verbose=1,
##D   control=list(max_treedepth=7),iter=150)
##D 
##D summary(f4,parmatrices=TRUE)
##D 
##D ctStanDiscretePars(f4,plot=TRUE) #auto and cross regressive plots over time
##D 
##D ctKalman(f4,timestep=.01,plot=TRUE,subjects=1,kalmanvec=c('y','etaprior'))
##D ctKalman(f4,timestep=.01,plot=TRUE,subjects=1:2,kalmanvec=c('y','etasmooth'))
##D ctKalman(f4,timestep=.01,plot=TRUE,subjects=1:2,kalmanvec=c('y','ysmooth'))
##D 
##D ctStanPlotPost(f4)
##D 
##D ctStanPostPredict(f4,wait=F)
##D 
##D 
##D 
##D #non-linear dedpendencies - based on m3 model (including intervention)
##D #specify intervention as dependent on extra parameter in PARS matrix, and latent process 1
##D 
##D m3nl <- ctModel(n.manifest = 1,n.latent = 1,n.TIpred = 1, type = 'stanct',
##D   manifestNames = c('Y1'), latentNames=c('L1'),TIpredNames = 'age',
##D   n.TDpred=1,TDpredNames = 'TD1', 
##D   TDPREDEFFECT=matrix(c('PARS[1,1] * state[1]'),nrow=1,ncol=1), 
##D   PARS=matrix(c('tdpredeffect'),1,1),
##D   DRIFT=matrix('drift11',nrow=1,ncol=1),
##D   DIFFUSION=matrix('diffusion11',nrow=1,ncol=1),
##D   CINT=matrix(c('cint1'),ncol=1),
##D   T0MEANS=matrix(c('t0m1'),ncol=1),
##D   T0VAR=matrix('t0var11',nrow=1,ncol=1),
##D   LAMBDA = diag(1),
##D   MANIFESTMEANS=matrix(0,ncol=1),
##D   MANIFESTVAR=matrix(c('merror1'),nrow=1,ncol=1))
##D 
##D m3nl$pars$indvarying[-which(m3nl$pars$matrix %in% 
##D   c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE
##D   
##D m3nl$pars$age_effect[-which(m3nl$pars$matrix %in% 
##D   c('T0MEANS','CINT','TDPREDEFFECT'))] <- FALSE
##D 
##D #here fit using optimization instead of sampling -- not appropriate in all cases!
##D f3nl <- ctStanFit(datalong = dat2, ctstanmodel = m3nl, 
##D   cores = setcores, chains = setchains,
##D   optimize=TRUE)
##D 
##D summary(f3nl)
##D 
##D #plot functions need updating for non-linearities! (as of ctsem v 2.7.3)
##D #extract can be used to extract samples and create own plots.
##D #last index of kalaman subobject denotes element of Kalman output.
##D # 1 = ll on std scale, 2= ll scale, 3=error, 4=prediction, 
##D # 5= eta prior, 6= eta upd
##D 
##D ctStanPostPredict(f3nl, datarows=1:100)
##D 
##D e=extract(f3nl)
##D subindex = which(f3nl$data$subject ==3) #specify subject
##D  
##D  matplot(f3nl$data$time[subindex], # Y predictions given earlier Y
##D    t(e$kalman[,subindex,4]), 
##D    type='l',lty=1,col=rgb(0,0,0,.1))
##D  
##D  points(f3nl$data$time[subindex], #actual Y
##D    f3nl$data$Y[subindex,1],type='p',col='red')
##D    
##D   matplot(f3nl$data$time[subindex], add = TRUE, #Generated Y from model
##D     t(e$Ygen[,subindex,1]), 
##D     type='l',lty=1,col=rgb(0,0,1,.05))
## End(Not run)



cleanEx()
nameEx("ctStanGenerateData")
### * ctStanGenerateData

flush(stderr()); flush(stdout())

### Name: ctStanGenerateData
### Title: Add a $generated object to ctstanfit object, with random data
###   generated from posterior of ctstanfit object
### Aliases: ctStanGenerateData

### ** Examples

## Not run: 
##D gen <- ctStanGenerateData(ctstantestfit, nsamples=3,fullposterior=TRUE)
##D plot(gen$generated$Y[,3,2],type='l') #Third random data sample, 2nd manifest var, all time points. 
## End(Not run)



cleanEx()
nameEx("ctStanKalman")
### * ctStanKalman

flush(stderr()); flush(stdout())

### Name: ctStanKalman
### Title: Get Kalman filter estimates from a ctStanFit object
### Aliases: ctStanKalman

### ** Examples

ctStanKalman(ctstantestfit)




cleanEx()
nameEx("ctStanModel")
### * ctStanModel

flush(stderr()); flush(stdout())

### Name: ctStanModel
### Title: Convert a frequentist (omx) ctsem model specification to
###   Bayesian (Stan).
### Aliases: ctStanModel

### ** Examples

model <- ctModel(type='omx', Tpoints=50,
n.latent=2, n.manifest=1, 
manifestNames='sunspots', 
latentNames=c('ss_level', 'ss_velocity'),
LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
DRIFT=matrix(c(0, 1,   'a21', 'a22'), nrow=2, ncol=2, byrow=TRUE),
MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
# MANIFESTVAR=matrix(0, nrow=1, ncol=1),
CINT=matrix(c(0, 0), nrow=2, ncol=1),
DIFFUSION=matrix(c(
  0, 0,
  0, "diffusion"), ncol=2, nrow=2, byrow=TRUE))

stanmodel=ctStanModel(model)





cleanEx()
nameEx("ctStanParMatrices")
### * ctStanParMatrices

flush(stderr()); flush(stdout())

### Name: ctStanParMatrices
### Title: Returns population system matrices from a ctStanFit object, and
###   vector of values for free parameters.
### Aliases: ctStanParMatrices

### ** Examples

## Not run: 
##D ctStanParMatrices(ctstantestfit,rnorm(17,0,.1))
## End(Not run)



cleanEx()
nameEx("ctStanParnames")
### * ctStanParnames

flush(stderr()); flush(stdout())

### Name: ctStanParnames
### Title: ctStanParnames
### Aliases: ctStanParnames

### ** Examples

ctStanParnames(ctstantestfit,substrings=c('pop_','popsd'))



cleanEx()
nameEx("ctStanPlotPost")
### * ctStanPlotPost

flush(stderr()); flush(stdout())

### Name: ctStanPlotPost
### Title: ctStanPlotPost
### Aliases: ctStanPlotPost

### ** Examples

ctStanPlotPost(ctstantestfit, rows=3:4)



cleanEx()
nameEx("ctStanPostPredict")
### * ctStanPostPredict

flush(stderr()); flush(stdout())

### Name: ctStanPostPredict
### Title: Compares model implied density and values to observed, for a
###   ctStanFit object.
### Aliases: ctStanPostPredict

### ** Examples

## Not run: 
##D ctStanPostPredict(ctstantestfit,wait=FALSE, shading=FALSE, datarows=1:25,diffsize=2)
## End(Not run)



cleanEx()
nameEx("ctStanTIpredMarginal")
### * ctStanTIpredMarginal

flush(stderr()); flush(stdout())

### Name: ctStanTIpredMarginal
### Title: Plot marginal relationships between covariates and parameters
###   for a ctStanFit object.
### Aliases: ctStanTIpredMarginal

### ** Examples

## Not run: 
##D ctStanTIpredMarginal(ctstantestfit,pars='CINT',tipred=3)
## End(Not run)



cleanEx()
nameEx("ctStanTIpredeffects")
### * ctStanTIpredeffects

flush(stderr()); flush(stdout())

### Name: ctStanTIpredeffects
### Title: Get time independent predictor effect estimates
### Aliases: ctStanTIpredeffects

### ** Examples

## Not run: 
##D #samples reduced here for speed
##D ctStanTIpredeffects(ctstantestfit,plot=TRUE,whichpars='CINT',nsamples=10,nsubjects=10)
## End(Not run)



cleanEx()
nameEx("ctStanUpdModel")
### * ctStanUpdModel

flush(stderr()); flush(stdout())

### Name: ctStanUpdModel
### Title: Update an already compiled and fit ctStanFit object
### Aliases: ctStanUpdModel

### ** Examples

## Not run: 
##D  newm<-ctModel(type='stanct',
##D   n.latent=ctstantestfit$ctstanmodel$n.latent,
##D   n.TDpred=ctstantestfit$ctstanmodel$n.TDpred,
##D   n.TIpred=ctstantestfit$ctstanmodel$n.TIpred,
##D   MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
##D   MANIFESTMEANS=matrix(0,nrow=ctstantestfit$ctstanmodel$n.manifest),
##D   CINT=matrix(c(0,'cint2'),ncol=1),
##D   n.manifest=ctstantestfit$ctstanmodel$n.manifest,
##D   LAMBDA=diag(2))
##D   
##D  newdat <- ctstantestdat
##D  newdat <- newdat[newdat[,'id']!=1,]
##D  newfit <- ctStanUpdModel(ctstantestfit, newdat, newm)
##D  
## End(Not run)



cleanEx()
nameEx("ctWideToLong")
### * ctWideToLong

flush(stderr()); flush(stdout())

### Name: ctWideToLong
### Title: ctWideToLong Convert ctsem wide to long format
### Aliases: ctWideToLong

### ** Examples

 #First load the example ctsem wide format data with absolute times
 data('datastructure')
 datastructure #contains two time intervals (dTx), therefore 3 time points.
 #Then convert to long format
 longexample <- ctWideToLong(datawide = datastructure, Tpoints=3, 
 n.manifest=3, manifestNames = c("Y1", "Y2", "Y3"),
 n.TDpred=1, TDpredNames = "TD1", 
 n.TIpred=2, TIpredNames = c("TI1", "TI2"))

 #Then convert the time intervals to absolute time
 long <- ctDeintervalise(datalong = longexample, id='id', dT='dT')
 long





cleanEx()
nameEx("ctstantestdat")
### * ctstantestdat

flush(stderr()); flush(stdout())

### Name: ctstantestdat
### Title: ctstantestdat
### Aliases: ctstantestdat

### ** Examples

## Not run: 
##D Tpoints=20
##D n.manifest=2
##D n.TDpred=1
##D n.TIpred=3
##D n.latent=2
##D n.subjects=5
##D gm<-ctModel(type='omx', Tpoints=Tpoints,n.latent=n.latent,
##D n.TDpred=n.TDpred,n.TIpred=n.TIpred,n.manifest=n.manifest,
##D   MANIFESTVAR=diag(0.5,2),
##D   TIPREDEFFECT=matrix(c(.5,0,0,-.7,0,2),nrow=2),
##D   TIPREDVAR=matrix(c(1,-.2,0, 0,1,0, 0,0,.5),nrow=3),
##D   TDPREDVAR=matrix(0,nrow=n.TDpred*(Tpoints),ncol=n.TDpred*(Tpoints)),
##D   TDPREDMEANS=matrix(round(exp(rnorm(n.TDpred*(Tpoints),-1.9,1)),0),
##D    nrow=n.TDpred*(Tpoints)),
##D    TDPREDEFFECT = matrix(c(1,-1),ncol=1),
##D   LAMBDA=diag(1,2),
##D   DRIFT=matrix(c(-.3,.2,0,-.2),nrow=2),
##D   DIFFUSION=matrix(c(2,1,0,2),2),
##D   CINT=matrix(c(0,0),nrow=2),
##D   T0MEANS=matrix(0,ncol=1,nrow=2),
##D   T0VAR=diag(100,2))
##D 
##D ctstantestdat<-ctGenerate(gm,n.subjects=n.subjects,burnin=30,
##D wide=FALSE,logdtsd=.4)
##D 
##D ctstantestdat[2,'Y1'] <- NA
##D ctstantestdat[ctstantestdat[,'id']==2,'TI1'] <- NA
##D ctstantestdat[2,'TD1'] <- NA
##D 
##D save(ctstantestdat,file='.\\data\\ctstantestdat.rda')
## End(Not run)



cleanEx()
nameEx("ctstantestfit")
### * ctstantestfit

flush(stderr()); flush(stdout())

### Name: ctstantestfit
### Title: ctstantestfit
### Aliases: ctstantestfit

### ** Examples

## Not run: 
##D #' 
##D checkm<-ctModel(type='stanct',
##D   n.latent=n.latent,n.TDpred=n.TDpred,n.TIpred=n.TIpred,
##D   MANIFESTVAR=matrix(c('merror',0,0,'merror'),2,2),
##D   MANIFESTMEANS=matrix(0,nrow=n.manifest),
##D   CINT=matrix(c('cint1','cint2'),ncol=1),
##D   n.manifest=n.manifest,LAMBDA=diag(2))
##D   
##D checkm$pars$indvarying[-c(7,13)] <- FALSE
##D checkm$pars$sdscale <- .1
##D  
##D checkm$pars[c(-1,-2, -21,-22) ,c('TI1_effect','TI2_effect','TI3_effect')] <- FALSE
##D 
##D ctstantestfit<-ctStanFit(ctstantestdat,checkm,iter=500, warmup=460,thin=2,chains=2,
##D   control=list(max_treedepth=5,adapt_delta=.8),save_warmup=FALSE)
##D summary(ctstantestfit)
##D save(ctstantestfit,file='.\\data\\ctstantestfit.rda')
##D paths <- sort(Sys.glob(c("data/*.rda", "data/*.RData")))
##D library(tools)
##D resaveRdaFiles(paths)
## End(Not run)



cleanEx()
nameEx("datastructure")
### * datastructure

flush(stderr()); flush(stdout())

### Name: datastructure
### Title: datastructure
### Aliases: datastructure

### ** Examples

## Not run: 
##D Tpoints=30
##D testm<-ctModel(Tpoints=Tpoints,n.latent=1,n.TDpred=1,n.TIpred=2,n.manifest=3,    
##D   LAMBDA=matrix(1,ncol=1,nrow=3),
##D   DRIFT=diag(-.3,1),
##D   DIFFUSION=diag(.1,1),
##D   CINT=diag(2,1),
##D   MANIFESTVAR=diag(1,3),
##D   TDPREDEFFECT=diag(.2,1),
##D   TIPREDEFFECT=matrix(.8,nrow=1,ncol=2),
##D   TDPREDVAR=diag(1,1*(Tpoints)),
##D   TIPREDVAR=diag(1,2)
##D )
##D longexample<-round(ctGenerate(testm,n.subjects=2,logdtsd = 1,burnin=3,wide=FALSE)[c(1:3,32:34),],2)
##D longexample[2,c(2,7)]<-NA
##D longexample[4,c(3)]<-NA
##D datastructure <- ctLongToWide(datalong = longexample,id='id',time='time',
##D   manifestNames = testm$manifestNames,TDpredNames = testm$TDpredNames,
##D   TIpredNames=testm$TIpredNames)
##D datastructure<-ctIntervalise(datawide = datastructure,
##D   Tpoints = 3,n.manifest = testm$n.manifest,n.TDpred = testm$n.TDpred,
##D   n.TIpred=testm$n.TIpred)
##D save(datastructure,file='.\\data\\datastructure.rda')
## End(Not run)



cleanEx()
nameEx("extract")
### * extract

flush(stderr()); flush(stdout())

### Name: extract
### Title: Extract samples from a ctStanFit object
### Aliases: extract

### ** Examples

e = extract(ctstantestfit)
head(e)



cleanEx()
nameEx("inv_logit")
### * inv_logit

flush(stderr()); flush(stdout())

### Name: inv_logit
### Title: Inverse logit
### Aliases: inv_logit

### ** Examples

inv_logit(-3)



cleanEx()
nameEx("isdiag")
### * isdiag

flush(stderr()); flush(stdout())

### Name: isdiag
### Title: Diagnostics for ctsem importance sampling
### Aliases: isdiag

### ** Examples

## Not run: 
##D #get data
##D sunspots<-sunspot.year
##D sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
##D id <- 1
##D time <- 1749:1924
##D datalong <- cbind(id, time, sunspots)
##D 
##D #setup model
##D model <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
##D  manifestNames='sunspots', 
##D  latentNames=c('ss_level', 'ss_velocity'),
##D   LAMBDA=matrix(c( -1, 'ma1' ), nrow=1, ncol=2),
##D   DRIFT=matrix(c(-.0001, 'a21', 1, 'a22'), nrow=2, ncol=2),
##D   MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
##D   CINT=matrix(c(0, 0), nrow=2, ncol=1),
##D   T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
##D   DIFFUSION=matrix(c(0.0001, 0, 0, "diffusion"), ncol=2, nrow=2))
##D 
##D model$pars$indvarying<-FALSE #Because single subject
##D model$pars$transform[14]<- '(param)*5+44 ' #Because not mean centered
##D model$pars$transform[4]<-'-log(exp(-param*1.5)+1)' #To avoid multi modality
##D 
##D #fit and plot importance sampling diagnostic
##D fit <- ctStanFit(datalong, model, chains=1,
##D   optimcontrol=list(isloops=5,finishsamples=500),optimize=TRUE)
##D isdiag(fit)
## End(Not run)



cleanEx()
nameEx("longexample")
### * longexample

flush(stderr()); flush(stdout())

### Name: longexample
### Title: longexample
### Aliases: longexample

### ** Examples

## Not run: 
##D #long example (using datastructure base)
##D Tpoints=30
##D testm<-ctModel(Tpoints=Tpoints,n.latent=1,n.TDpred=1,n.TIpred=2,n.manifest=3,    
##D   LAMBDA=matrix(1,ncol=1,nrow=3),
##D   DRIFT=diag(-.3,1),
##D   DIFFUSION=diag(.1,1),
##D   CINT=diag(2,1),
##D   MANIFESTVAR=diag(1,3),
##D   TDPREDEFFECT=diag(.2,1),
##D   TIPREDEFFECT=matrix(.8,nrow=1,ncol=2),
##D   TDPREDVAR=diag(1,1*(Tpoints)),
##D   TIPREDVAR=diag(1,2)
##D )
##D longexample<-round(ctGenerate(testm,n.subjects=2,logdtsd = 1,burnin=3,wide=FALSE)[c(1:3,32:35),],2)
##D longexample[2,c(2,7)]<-NA
##D longexample[4,c(3)]<-NA
##D longexample
##D save(longexample,file='.\\data\\longexample.rda')
## End(Not run)



cleanEx()
nameEx("msquare")
### * msquare

flush(stderr()); flush(stdout())

### Name: msquare
### Title: Right multiply a matrix by its transpose.
### Aliases: msquare

### ** Examples

msquare(t(chol(diag(3,4)+1)))



cleanEx()
nameEx("optimstan")
### * optimstan

flush(stderr()); flush(stdout())

### Name: optimstan
### Title: Optimize / importance sample a stan or ctStan model.
### Aliases: optimstan

### ** Examples

 sunspots<-sunspot.year
 sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
 id <- 1
 time <- 1749:1924
datalong <- cbind(id, time, sunspots)

#setup model
 ssmodel <- ctModel(type='stanct', n.latent=2, n.manifest=1,
  manifestNames='sunspots',
  latentNames=c('ss_level', 'ss_velocity'),
   LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
   DRIFT=matrix(c(0, 'a21', 1, 'a22'), nrow=2, ncol=2),
   MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
   MANIFESTVAR=diag(0,1),
   CINT=matrix(c(0, 0), nrow=2, ncol=1),
   T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
   DIFFUSION=matrix(c(0, 0, 0, "diffusion"), ncol=2, nrow=2))

 ssmodel$pars$indvarying<-FALSE #Because single subject
 ssmodel$pars$offset[14]<- 44 #Because not mean centered
 ssmodel$pars[4,c('transform','offset')]<- c(1,0) #To avoid multi modality

#fit using optimization without importance sampling
ssfit <- ctStanFit(datalong[1:50,], #limited data for example
  ssmodel, optimize=TRUE,optimcontrol=list(deoptim=FALSE,isloops=0,finishsamples=50))

#output
summary(ssfit)



cleanEx()
nameEx("plot.ctStanFit")
### * plot.ctStanFit

flush(stderr()); flush(stdout())

### Name: plot.ctStanFit
### Title: plot.ctStanFit
### Aliases: plot.ctStanFit ctStanPlot

### ** Examples

## Not run: 
##D plot(ctstantestfit,types=c('regression','kalman','priorcheck'))
##D 
##D ### complete example
##D plot(ctstantestfit)
##D 
##D #### example plot using rstan functions
##D rstan::stan_trace(ctstantestfit$stanfit, 
##D pars=ctStanParnames(ctstantestfit,'pop_DRIFT'))
## End(Not run)



cleanEx()
nameEx("plot.ctStanModel")
### * plot.ctStanModel

flush(stderr()); flush(stdout())

### Name: plot.ctStanModel
### Title: Prior plotting
### Aliases: plot.ctStanModel

### ** Examples

model <- ctModel(type='omx', Tpoints=50,
n.latent=2, n.manifest=1, 
manifestNames='sunspots', 
latentNames=c('ss_level', 'ss_velocity'),
LAMBDA=matrix(c( 1, 'ma1' ), nrow=1, ncol=2),
DRIFT=matrix(c(0, 1,   'a21', 'a22'), nrow=2, ncol=2, byrow=TRUE),
MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
# MANIFESTVAR=matrix(0, nrow=1, ncol=1),
CINT=matrix(c(0, 0), nrow=2, ncol=1),
DIFFUSION=matrix(c(
  0, 0,
  0, "diffusion"), ncol=2, nrow=2, byrow=TRUE))

stanmodel=ctStanModel(model)
plot(stanmodel,rows=8)



cleanEx()
nameEx("plot.ctsemFit")
### * plot.ctsemFit

flush(stderr()); flush(stdout())

### Name: plot.ctsemFit
### Title: Plotting function for object class ctsemFit
### Aliases: plot.ctsemFit

### ** Examples

## Examples set to 'dontrun' because they take longer than 5s.

### example from Driver, Oud, Voelkle (2015), 
### simulated happiness and leisure time with unobserved heterogeneity.
## Not run: 
##D data(ctExample1)
##D traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
##D   manifestNames=c('LeisureTime', 'Happiness'), 
##D   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
##D traitfit <- ctFit(datawide=ctExample1, ctmodelobj=traitmodel)
##D plot(traitfit, wait=FALSE)
## End(Not run)



cleanEx()
nameEx("plot.ctsemFitMeasure")
### * plot.ctsemFitMeasure

flush(stderr()); flush(stdout())

### Name: plot.ctsemFitMeasure
### Title: Misspecification plot using ctCheckFit output
### Aliases: plot.ctsemFitMeasure

### ** Examples

## Not run: 
##D data(ctExample1)
##D traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
##D   manifestNames=c('LeisureTime', 'Happiness'), 
##D   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
##D traitfit <- ctFit(datawide=ctExample1, ctmodelobj=traitmodel)
##D 
##D check <- ctCheckFit(traitfit,niter=5)
##D plot(check)
## End(Not run)



cleanEx()
nameEx("sdpcor2cov")
### * sdpcor2cov

flush(stderr()); flush(stdout())

### Name: sdpcor2cov
### Title: sdcor2cov
### Aliases: sdpcor2cov

### ** Examples

testmat <- diag(exp(rnorm(5,-3,2)),5) #generate arbitrary std deviations
testmat[row(testmat) > col(testmat)] <- runif((5^2-5)/2, -1, 1) 
print(testmat)
covmat <- sdpcor2cov(testmat) #convert to covariance
cov2cor(covmat) #convert covariance to correlation



cleanEx()
nameEx("stanWplot")
### * stanWplot

flush(stderr()); flush(stdout())

### Name: stanWplot
### Title: Runs stan, and plots sampling information while sampling.
### Aliases: stanWplot

### ** Examples

## Not run: 
##D #### example 1 
##D scode <- "
##D parameters {
##D   real y[2]; 
##D } 
##D model {
##D   y[1] ~ normal(0, .5);
##D   y[2] ~ double_exponential(0, 2);
##D } 
##D "
##D sm <- stan_model(model_code = scode)
##D fit1 <- stanWplot(object = sm,iter = 100000,chains=2,cores=1)
## End(Not run)



cleanEx()
nameEx("stan_checkdivergences")
### * stan_checkdivergences

flush(stderr()); flush(stdout())

### Name: stan_checkdivergences
### Title: Analyse divergences in a stanfit object
### Aliases: stan_checkdivergences

### ** Examples

## Not run: 
##D stan_checkdivergences(myfitobj)
## End(Not run)



cleanEx()
nameEx("stan_confidenceRegion")
### * stan_confidenceRegion

flush(stderr()); flush(stdout())

### Name: stan_confidenceRegion
### Title: Extract functions of multiple variables from a stanfit object
### Aliases: stan_confidenceRegion

### ** Examples

temp<-stan_confidenceRegion(stanfit=ctstantestfit$stanfit, 
  parstrings=c('pop_DRIFT[1,2]','pop_DRIFT[2,1]'))
t(apply(temp,2,quantile))



cleanEx()
nameEx("stan_postcalc")
### * stan_postcalc

flush(stderr()); flush(stdout())

### Name: stan_postcalc
### Title: Compute functions of matrices from samples of a stanfit object
### Aliases: stan_postcalc

### ** Examples

temp<-stan_postcalc(stanfit=ctstantestfit$stanfit, 
  object='DRIFT', objectindices='all', calc='exp(object)')



cleanEx()
nameEx("stan_unconstrainsamples")
### * stan_unconstrainsamples

flush(stderr()); flush(stdout())

### Name: stan_unconstrainsamples
### Title: Convert samples from a stanfit object to the unconstrained scale
### Aliases: stan_unconstrainsamples

### ** Examples

## Not run: 
##D umat <- stan_unconstrainsamples(ctstantestfit$stanfit, ctstantestfit$standata)
## End(Not run)



cleanEx()
nameEx("summary.ctStanFit")
### * summary.ctStanFit

flush(stderr()); flush(stdout())

### Name: summary.ctStanFit
### Title: summary.ctStanFit
### Aliases: summary.ctStanFit

### ** Examples

summary(ctstantestfit)



cleanEx()
nameEx("summary.ctsemFit")
### * summary.ctsemFit

flush(stderr()); flush(stdout())

### Name: summary.ctsemFit
### Title: Summary function for ctsemFit object
### Aliases: summary.ctsemFit

### ** Examples

## Examples set to 'dontrun' because they take longer than 5s. 

## Not run: 
##D ### example from Driver, Oud, Voelkle (2015), 
##D ### simulated happiness and leisure time with unobserved heterogeneity.
##D data(ctExample1)
##D traitmodel <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2), 
##D   manifestNames=c('LeisureTime', 'Happiness'), 
##D   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
##D traitfit <- ctFit(datawide=ctExample1, ctmodelobj=traitmodel)
##D summary(traitfit,timeInterval=1)
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
