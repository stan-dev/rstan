pkgname <- "clinDR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('clinDR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Brextract")
### * Brextract

flush(stderr()); flush(stdout())

### Name: "Extract.emaxsim"
### Title: Extract a simulation from the output of emaxsim
### Aliases: "Extract.emaxsim" [.emaxsim
### Keywords: nonlinear

### ** Examples


## Not run: 
##D ## code change random number seed
##D 
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D D1 <- emaxsim(nsim,gen.parm,modType=3)
##D e49<-D1[49]                  #### extract 49th simulation
##D 
## End(Not run)
## Don't show: 
## code change random number seed

doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  

D1 <- emaxsim(nsim=2,gen.parm,modType=3,nproc=1)
e49<-D1[2]                  
## End(Don't show)



cleanEx()
nameEx("BrextractB")
### * BrextractB

flush(stderr()); flush(stdout())

### Name: "Extract.emaxsimB"
### Title: Extract a simulation from the output of emaxsimB
### Aliases: "Extract.emaxsimB" [.emaxsimB
### Keywords: nonlinear

### ** Examples


## Not run: 
##D 
##D save.seed<-.Random.seed
##D set.seed(12357)
##D 
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D Ndose<-length(doselev)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-4.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
##D 		sigmaup=30,edDF=5)
##D mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)
##D 
##D D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)
##D 
##D out<-D1[2]
##D 
##D 
##D .Random.seed<-save.seed
## End(Not run)



cleanEx()
nameEx("DRDensityPlot")
### * DRDensityPlot

flush(stderr()); flush(stdout())

### Name: DRDensityPlot
### Title: Plot Bayes or confidence interval density contours over a grid
###   of points (usually dose or time)
### Aliases: DRDensityPlot
### Keywords: nonlinear hplot

### ** Examples


## Not run: 
##D data(examples14)
##D exdat<-examples14[[32]]
##D 
##D fitout<-fitEmax(exdat$y,exdat$dose,modType=3,count=exdat$nsize,
##D 								msSat=(exdat$sd)^2)
##D 
##D dgrid<-seq(0,1,length=100)
##D seout95<-predict(fitout,dgrid,clev=0.95)
##D seout90<-predict(fitout,dgrid,clev=0.9)
##D seout80<-predict(fitout,dgrid,clev=0.8)
##D seout50<-predict(fitout,dgrid,clev=0.5)
##D 
##D qlev<-c(0.025,0.05,0.10,0.25)
##D 
##D qL<-cbind(seout95$ubdif,seout90$ubdif,seout80$ubdif,seout50$ubdif)
##D qH<-cbind(seout95$lbdif,seout90$lbdif,seout80$lbdif,seout50$lbdif)
##D 
##D DRDensityPlot(dgrid,qL,qH,qlevL=qlev,xlab='Dose',ylab='Diff with PBO')
##D 
## End(Not run)
## Don't show: 
data(examples14)
exdat<-examples14[[32]]

dgrid<-seq(0,1,length=5)

qlev<-c(0.10,0.25)

qL<-matrix(c(0.000000, 0.000000,
1.181590, 1.093189,
1.301505, 1.220726,
1.354046, 1.273955,
1.384266, 1.303586),ncol=2,byrow=TRUE)

qH<-matrix(c(0.0000000, 0.0000000,
0.8083449, 0.8967468,
0.9604440, 1.0412232,
1.0158898, 1.0959808,
1.0436238, 1.1243036),ncol=2,byrow=TRUE)

DRDensityPlot(dgrid,qL,qH,qlevL=qlev,xlab='Dose',ylab='Diff with PBO')
## End(Don't show)



cleanEx()
nameEx("FixedMean")
### * FixedMean

flush(stderr()); flush(stdout())

### Name: FixedMean
### Title: Fixed means (proportions) random data constructor for emaxsim
###   for continuous or binary data
### Aliases: FixedMean
### Keywords: nonlinear

### ** Examples


## Not run: 
##D ##  example changes the random number seed
##D 
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)  
##D   
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen.parm<-FixedMean(n,doselev,meanlev,sdy,pop)  
##D 
##D ### 4-parameter example
##D 
##D est<-c( log(6.67636), 2,0, -3.18230)
##D doselev<-c(0,5,10,25,50,150)
##D meanlev<-emaxfun(doselev,est)
##D 
##D gen.parm4<-FixedMean(n=c(99,95,98,94,98,99),doselev,
##D                      meanlev,resSD=3.87,parm=est)
##D D4 <- emaxsim(nsim=100,gen.parm4,modType=4,negEmax=TRUE)
##D summary(D4)
## End(Not run)
## Don't show: 
##  example changes the random number seed

doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)  
  
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy,pop)  

## End(Don't show)



cleanEx()
nameEx("RandEmax")
### * RandEmax

flush(stderr()); flush(stdout())

### Name: RandEmax
### Title: Random data constructor function for emaxsim creating random
###   parameters for an Emax model for continuous or binary data.
### Aliases: RandEmax
### Keywords: nonlinear

### ** Examples

simParm<-RandEmax(n=c(99,95,98,94,98,98),doselev=c(0,5,10,25,50,150),
         parmE0=c(-2.6,2.5),p50=25,parmEmax=c(-1.25,2),resSD=3.88)



cleanEx()
nameEx("SeEmax")
### * SeEmax

flush(stderr()); flush(stdout())

### Name: SeEmax
### Title: Asymptotic SE for dose response estimates from a 3- or 4-
###   parameter Emax model
### Aliases: SeEmax
### Keywords: nonlinear

### ** Examples


## Not run: 
##D 
##D ## this example changes the random number seed
##D doselev<-c(0,5,25,50,100,250)
##D n<-c(78,81,81,81,77,80)
##D dose<-rep(doselev,n)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D led50<-log(ed50)
##D emax<-15.127726
##D lambda=1.8
##D sdy<-7.967897
##D pop<-c(led50=led50,lambda=lambda,emax=emax,e0=e0)    
##D meanresp<-emaxfun(dose,pop)  
##D y<-rnorm(sum(n),meanresp,sdy)
##D nls.fit<-nls(y ~ e0 + (emax * dose^lambda)/(dose^lambda + exp(led50*lambda)), 
##D                          start = pop, control = nls.control(
##D                          maxiter = 100),trace=TRUE,na.action=na.omit)
##D 
##D 
##D SeEmax(nls.fit,doselev=c(60,120),modType=4)
##D SeEmax(list(coef(nls.fit),vcov(nls.fit)),c(60,120),modType=4)
## End(Not run)
## Don't show: 

## this example changes the random number seed
doselev<-c(0,5,25,50,100,250)
n<-c(78,81,81,81,77,80)
dose<-rep(doselev,n)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
led50<-log(ed50)
emax<-15.127726
lambda=1.8
sdy<-7.967897
pop<-c(led50=led50,lambda=lambda,emax=emax,e0=e0)    
meanresp<-emaxfun(dose,pop)  
y<-rnorm(sum(n),meanresp,sdy)
nls.fit<-nls(y ~ e0 + (emax * dose^lambda)/(dose^lambda + exp(led50*lambda)), 
                         start = pop, control = nls.control(
                         maxiter = 100),trace=TRUE,na.action=na.omit)


SeEmax(nls.fit,doselev=c(60,120),modType=4)
## End(Don't show)



cleanEx()
nameEx("checkMonoEmax")
### * checkMonoEmax

flush(stderr()); flush(stdout())

### Name: checkMonoEmax
### Title: Bayes posterior predictive test for Emax (monotone) model fit
### Aliases: checkMonoEmax
### Keywords: nonlinear

### ** Examples

## Not run: 
##D 
##D data("examples14")
##D exdat<-examples14[[6]]
##D 
##D prior<-prior.control(epmu=0,epsd=10,emaxmu=0,emaxsd=10,p50=0.25,
##D 				sigmalow=0.01,sigmaup=3)
##D mcmc<-mcmc.control(chains=3)
##D 
##D fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
##D 				count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)
##D parms<-coef(fitout)[,1:4]  #use first intercept
##D 
##D checkMonoEmax(y=exdat$y, dose=exdat$dose, parm=parms, sigma2=(sigma(fitout))^2,
##D       nvec=exdat$nsize, trend='negative')
##D       
## End(Not run)
## Don't show: 

data("examples14")
exdat<-examples14[[6]]


parms<-cbind(rnorm(5,-.69,.2),rnorm(5,.5,.01),rnorm(5,-5.19,0.2),rnorm(5,1.8,.02))
sig2<-rnorm(5,0.4,0.01)

checkMonoEmax(y=exdat$y, dose=exdat$dose, parm=parms, sigma2=sig2^2,
      nvec=exdat$nsize, trend='negative')
## End(Don't show)



cleanEx()
nameEx("coef")
### * coef

flush(stderr()); flush(stdout())

### Name: coefEmax
### Title: Extract Emax model parameter estimates
### Aliases: coef.fitEmax coef.fitEmaxB coef.emaxsim coef.emaxsimB
### Keywords: nonlinear

### ** Examples

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)
coef(testout)



cleanEx()
nameEx("emaxalt")
### * emaxalt

flush(stderr()); flush(stdout())

### Name: emaxalt
### Title: Fit 4- or 3-parameter Emax model substituting simpler curves if
###   convergence not achieved.
### Aliases: emaxalt
### Keywords: nonlinear

### ** Examples


save.seed<-.Random.seed
set.seed(12357)

doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
dose<-rep(doselev,n)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanresp<-emaxfun(dose,pop)  
y<-rnorm(sum(n),meanresp,sdy)

simout<-emaxalt(y,dose)

simout2<-emaxalt(y,dose,modType=4)

.Random.seed<-save.seed



cleanEx()
nameEx("emaxfun")
### * emaxfun

flush(stderr()); flush(stdout())

### Name: emaxfun
### Title: Vectorized versions of the hyperbolic and sigmoidal Emax models
### Aliases: emaxfun
### Keywords: nonlinear

### ** Examples


doselev<-c(0,5,25,50,100)
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
lambda=2
parm<-c(log(ed50),lambda,emax,e0)
plot(doselev,emaxfun(doselev,parm))




cleanEx()
nameEx("emaxsim")
### * emaxsim

flush(stderr()); flush(stdout())

### Name: emaxsim
### Title: Simulate Emax maximum likelihood estimation
### Aliases: emaxsim
### Keywords: nonlinear

### ** Examples


## Not run: 
##D ## emaxsim changes the random number seed
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D Ndose<-length(doselev)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-4.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D D1 <- emaxsim(nsim,gen,modType=3)
##D summary(D1,testalph=0.05)
##D 
##D D4 <- emaxsim(nsim,gen,modType=4)
##D summary(D4,testalph=0.05)
## End(Not run)
## Don't show: 
## emaxsim changes the random number seed
nsim<-2
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
Ndose<-length(doselev)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-4.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen<-FixedMean(n,doselev,meanlev,sdy)  

D1 <- emaxsim(nsim,gen,modType=3,nproc=1)
## End(Don't show)



cleanEx()
nameEx("emaxsimB")
### * emaxsimB

flush(stderr()); flush(stdout())

### Name: emaxsimB
### Title: Simulate Emax Bayesian estimation
### Aliases: emaxsimB
### Keywords: Bayes Emax

### ** Examples

## Not run: 
##D 
##D ### emaxsimB changes the random number seed
##D 
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D Ndose<-length(doselev)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-4.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
##D 		sigmaup=30,edDF=5)
##D mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,
##D 		propInit=0.15,adapt_delta = 0.95)
##D 
##D ### custom code to compute the distribution of the dose yielding
##D ### a target diff with pbo
##D customCode<-function(parms,residSD,pVal,dose,y,customParms){
##D 	target<-customParms
##D 	ed50<-exp(parms[,1])
##D 	emax<-parms[,2]
##D 	td<-ifelse(emax-target>0,ed50*(target/(emax-target)),Inf)
##D 	tdest<-median(td)
##D 	lb<-quantile(td,0.1)
##D 	ub<-quantile(td,0.9)
##D 	return(c(td=tdest,lb=lb,ub=ub))
##D }
##D 
##D D1 <- emaxsimB(nsim,gen, prior, modType=3,seed=12357,mcmc=mcmc,check=FALSE,
##D 				customCode=customCode,customParms=1.0)
##D D1
## End(Not run)



cleanEx()
nameEx("examples14")
### * examples14

flush(stderr()); flush(stdout())

### Name: examples14
### Title: Data from 2014 dose response meta-analysis
### Aliases: examples14 meta14
### Keywords: datasets

### ** Examples

data(examples14)
objects(examples14[[1]])
names(meta14)



cleanEx()
nameEx("examples16")
### * examples16

flush(stderr()); flush(stdout())

### Name: examples16
### Title: Data from 2016 dose response meta-analysis
### Aliases: examples16 meta16
### Keywords: datasets

### ** Examples

data(examples16)
objects(examples16[[1]])
names(meta16)



cleanEx()
nameEx("fitEmax")
### * fitEmax

flush(stderr()); flush(stdout())

### Name: fitEmax
### Title: ML fit of hyperbolic or sigmoidal Emax models to
###   continuous/binary dose response data.
### Aliases: fitEmax
### Keywords: nonlinear

### ** Examples


## the example changes the random number seed

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)



cleanEx()
nameEx("fitEmaxB")
### * fitEmaxB

flush(stderr()); flush(stdout())

### Name: fitEmaxB
### Title: Bayesian fit of hyperbolic or sigmoidal Emax models to
###   continuous/binary dose response data.
### Aliases: fitEmaxB
### Keywords: dose response mcmc

### ** Examples

## Not run: 
##D 
##D data("examples14")
##D exdat<-examples14[[1]]
##D 
##D prior<-prior.control(epmu=0,epsd=4,emaxmu=0,emaxsd=4,p50=0.1,
##D 				sigmalow=0.01,sigmaup=3)
##D 										
##D mcmc<-mcmc.control(chains=3)
##D 
##D fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
##D 				count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)
## End(Not run)



cleanEx()
nameEx("nllogis")
### * nllogis

flush(stderr()); flush(stdout())

### Name: nllogis
### Title: The negative log likelihood function for a 3- or 4- parameter
###   Emax model on the logit scale for binary dose response.
### Aliases: nllogis
### Keywords: nonlinear

### ** Examples

data(examples14)
with(examples14[[8]],nllogis(parms=c(log(.17),-3.26,-0.15), y, dose))



cleanEx()
nameEx("plot.emaxsim")
### * plot.emaxsim

flush(stderr()); flush(stdout())

### Name: plot.emaxsim
### Title: Plot the output of emaxsim
### Aliases: plot.emaxsim
### Keywords: nonlinear

### ** Examples


## Not run: 
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop.parm<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop.parm)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
##D D1 <- emaxsim(nsim,gen.parm)
##D 
##D plot(D1,id=3)
## End(Not run)
## Don't show: 
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop.parm<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop.parm)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim=3,gen.parm,nproc=1)

plot(D1,id=3)
## End(Don't show)



cleanEx()
nameEx("plot.emaxsimB")
### * plot.emaxsimB

flush(stderr()); flush(stdout())

### Name: plot.emaxsimB
### Title: Plot the output of emaxsimB
### Aliases: plot.emaxsimB
### Keywords: Bayes Emax

### ** Examples

## Not run: 
##D ## emaxsimB changes the random number seeds
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D Ndose<-length(doselev)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-4.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
##D 		sigmaup=30,edDF=5)
##D mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)
##D 
##D D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)
##D 
##D plot(D1,id=3)
## End(Not run)



cleanEx()
nameEx("plot.emaxsimBobj")
### * plot.emaxsimBobj

flush(stderr()); flush(stdout())

### Name: plot.emaxsimBobj
### Title: Plot dose response from a data set generated by emaxsimB
### Aliases: plot.emaxsimBobj

### ** Examples

## Not run: 
##D 
##D ## emaxsimB changes the random number seed
##D 
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D Ndose<-length(doselev)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-4.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D 
##D prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
##D 		sigmaup=30,edDF=5)
##D mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)
##D 
##D D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)
##D 
##D plot(D1[2])
## End(Not run)



cleanEx()
nameEx("plot.emaxsimobj")
### * plot.emaxsimobj

flush(stderr()); flush(stdout())

### Name: plot.emaxsimobj
### Title: Plot dose response from a data set generated by emaxsim
### Aliases: plot.emaxsimobj
### Keywords: nonlinear

### ** Examples

## Not run: 
##D ## emaxsim changes the random number seed
##D 
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
##D D1 <- emaxsim(nsim,gen.parm)
##D e49<-D1[49]
##D 
##D plot(e49,clev=0.8)  
## End(Not run)
## Don't show: 
## emaxsim changes the random number seed

nsim<-2
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen.parm,nproc=1)
e2<-D1[2]

plot(e2,clev=0.8)  
## End(Don't show)



cleanEx()
nameEx("plot.fitEmax")
### * plot.fitEmax

flush(stderr()); flush(stdout())

### Name: plot.fitEmax
### Title: Plot a Emax model and dose group means.
### Aliases: plot.fitEmax
### Keywords: nonlinear

### ** Examples

### example changes the random number seed

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop.parm<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop.parm)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)

plot(testout)



cleanEx()
nameEx("plot.fitEmaxB")
### * plot.fitEmaxB

flush(stderr()); flush(stdout())

### Name: plot.fitEmaxB
### Title: Plot a Emax model and dose group means.
### Aliases: plot.fitEmaxB
### Keywords: nonlinear

### ** Examples

## Not run: 
##D data("examples14")
##D exdat<-examples14[[1]]
##D 
##D prior<-prior.control(epmu=0,epsd=4,emaxmu=0,emaxsd=4,p50=0.1,
##D 				sigmalow=0.01,sigmaup=3)
##D mcmc<-mcmc.control(chains=3)
##D 
##D fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
##D 				count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)
##D 				
##D plot(fitout)
## End(Not run)



cleanEx()
nameEx("plot.plotB")
### * plot.plotB

flush(stderr()); flush(stdout())

### Name: plot.plotB
### Title: Plot Bayes dose response curve and dose group means
### Aliases: plot.plotB
### Keywords: nonlinear hplot

### ** Examples

## Not run: 
##D data("examples14")
##D exdat<-examples14[[6]]
##D 
##D prior<-prior.control(epmu=0,epsd=10,emaxmu=0,emaxsd=10,p50=0.25,
##D 				sigmalow=0.01,sigmaup=3)
##D mcmc<-mcmc.control(chains=3)
##D 
##D fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
##D 				count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)
##D parms<-coef(fitout)[1:4]   ### use first intercept
##D 
##D outB<-plotB(exdat$y,exdat$dose,parms, sigma2=(sigma(fitout))^2,
##D 			ylab="Change in EDD")
##D 			
##D plot(outB,plotDif=TRUE)
## End(Not run)
## Don't show: 
data("examples14")
exdat<-examples14[[6]]

parms<-matrix(c(-0.1665350, 0.3657811, -5.660137, 1.744753, 0.4050860,
-0.8463137, 0.3837361, -4.877676, 1.784098, 0.3943782,
-1.1811274, 0.3767222, -4.921861, 1.873861, 0.4266011,
 0.4729616, 0.3157714, -6.322768, 1.780517, 0.3646588,
 0.4255880, 0.3336959, -6.251558, 1.775438, 0.3657461),ncol=5,byrow=TRUE)

outB<-plotB(exdat$y,exdat$dose,parms[,1:4], sigma2=(parms[,5])^2,
			ylab="Change in EDD")
			
plot(outB,plotDif=TRUE)
## End(Don't show)



cleanEx()
nameEx("plotB")
### * plotB

flush(stderr()); flush(stdout())

### Name: plotB
### Title: Plot Bayes dose response curve and dose group means
### Aliases: plotB
### Keywords: nonlinear hplot

### ** Examples

## Not run: 
##D data("examples14")
##D exdat<-examples14[[6]]
##D 
##D prior<-prior.control(epmu=0,epsd=10,emaxmu=0,emaxsd=10,p50=0.25,
##D 				sigmalow=0.01,sigmaup=3)
##D mcmc<-mcmc.control(chains=3)
##D 
##D fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
##D 				count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)
##D parms<-coef(fitout)[1:4]  ## use first intercept
##D 
##D outB<-plotB(exdat$y,exdat$dose,parms, sigma2=(sigma(fitout))^2,
##D 			ylab="Change in EDD")
##D 			
##D plot(outB,plotDif=TRUE)
## End(Not run)
## Don't show: 
data("examples14")
exdat<-examples14[[6]]

parms<-matrix(c(-0.1665350, 0.3657811, -5.660137, 1.744753, 0.4050860,
-0.8463137, 0.3837361, -4.877676, 1.784098, 0.3943782,
-1.1811274, 0.3767222, -4.921861, 1.873861, 0.4266011,
 0.4729616, 0.3157714, -6.322768, 1.780517, 0.3646588,
 0.4255880, 0.3336959, -6.251558, 1.775438, 0.3657461),ncol=5,byrow=TRUE)

outB<-plotB(exdat$y,exdat$dose,parms[,1:4], sigma2=(parms[,5])^2,
			ylab="Change in EDD")
			
plot(outB,plotDif=TRUE)
## End(Don't show)



cleanEx()
nameEx("plotBdensity")
### * plotBdensity

flush(stderr()); flush(stdout())

### Name: plotBdensity
### Title: Density plot displaying Bayes prior or posterior dose response
### Aliases: plotBdensity
### Keywords: nonlinear hplot

### ** Examples


## Not run: 
##D data("examples14")
##D exdat<-examples14[[6]]
##D 
##D prior<-prior.control(epmu=0,epsd=10,emaxmu=0,emaxsd=10,p50=0.25,
##D 				sigmalow=0.01,sigmaup=3)
##D mcmc<-mcmc.control(chains=3)
##D 
##D fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
##D 				count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)
##D parms<-coef(fitout)[1:4]   ## use first intercept
##D 
##D dgrid<-seq(0,1,length=100)
##D 
##D plotBdensity(dgrid,parm=parms)
##D 
##D plotBdensity(dgrid,parm=parms,plotDif=TRUE,
##D        xlab='Dose',ylab='Dif with PBO')
## End(Not run)
## Don't show: 
data("examples14")
exdat<-examples14[[6]]

parms<-matrix(c(-0.1665350, 0.3657811, -5.660137, 1.744753, 0.4050860,
-0.8463137, 0.3837361, -4.877676, 1.784098, 0.3943782,
-1.1811274, 0.3767222, -4.921861, 1.873861, 0.4266011,
 0.4729616, 0.3157714, -6.322768, 1.780517, 0.3646588,
 0.4255880, 0.3336959, -6.251558, 1.775438, 0.3657461),ncol=5,byrow=TRUE)

dgrid<-seq(0,1,length=5)

plotBdensity(dgrid,parm=parms[,1:4])
## End(Don't show)



cleanEx()
nameEx("plotD")
### * plotD

flush(stderr()); flush(stdout())

### Name: plotD
### Title: Basic plot of dose group means
### Aliases: plotD
### Keywords: hplot

### ** Examples

data(examples14)
with(examples14[[2]],plotD(y,dose,meansOnly=TRUE,se=TRUE,sem=sem,ylab=
"Y",xlab="Dose(mg)"))



cleanEx()
nameEx("predict.emaxalt")
### * predict.emaxalt

flush(stderr()); flush(stdout())

### Name: predict.emaxalt
### Title: Mean response and SE for specified doses for a simulated object
###   output by function emaxalt
### Aliases: predict.emaxalt
### Keywords: nonlinear

### ** Examples


## Not run: 
##D ## random number seed changed by this example
##D 
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D dose<-rep(doselev,n)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop.parm<-c(log(ed50),e0,emax)    
##D meanresp<-emaxfun(dose,pop.parm)  
##D y<-rnorm(sum(n),meanresp,sdy)
##D 
##D simout<-emaxalt(y,dose)
##D predict(simout,c(75,150))
##D 
##D simout2<-emaxalt(y,dose,modType=4)
##D predict(simout2,c(75,150))
## End(Not run)
## Don't show: 
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)
dose<-rep(doselev,n)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop.parm<-c(log(ed50),e0,emax)    
meanresp<-emaxfun(dose,pop.parm)  
y<-rnorm(sum(n),meanresp,sdy)

simout<-emaxalt(y,dose)
predict(simout,c(75,150))
## End(Don't show)



cleanEx()
nameEx("predict.emaxsim")
### * predict.emaxsim

flush(stderr()); flush(stdout())

### Name: predict.emaxsim
### Title: Mean response and SE for specified doses for each replicate data
###   set in an emaxsim object
### Aliases: predict.emaxsim
### Keywords: nonlinear

### ** Examples


## Not run: 
##D ## random number seed changed by this example
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop.parm<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop.parm)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
##D D1 <- emaxsim(nsim,gen.parm)
##D 
##D predout<-predict(D1,c(75,150))
## End(Not run)
## Don't show: 
## random number seed changed by this example
nsim<-3
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop.parm<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop.parm)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen.parm,nproc=1)

predout<-predict(D1,c(75,150))
## End(Don't show)



cleanEx()
nameEx("predict.emaxsimB")
### * predict.emaxsimB

flush(stderr()); flush(stdout())

### Name: predict.emaxsimB
### Title: Mean response and SE for each replicate data set in an emaxsimB
###   object
### Aliases: predict.emaxsimB
### Keywords: nonlinear

### ** Examples


## Not run: 
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D Ndose<-length(doselev)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-4.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
##D 		sigmaup=30,edDF=5)
##D mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,
##D 		propInit=0.15,adapt_delta = 0.95)
##D 
##D 
##D D1 <- emaxsimB(nsim,gen, prior, modType=3,seed=12357,mcmc=mcmc,check=FALSE)
##D 
##D predict(D1,dose=20)
## End(Not run)



cleanEx()
nameEx("predict.emaxsimBobj")
### * predict.emaxsimBobj

flush(stderr()); flush(stdout())

### Name: predict.emaxsimBobj
### Title: Mean response estimates (posterior means) and SE (posterior SD)
###   for specified doses for a simulated emaxsimBobj object
### Aliases: predict.emaxsimBobj
### Keywords: nonlinear

### ** Examples


## Not run: 
##D ### emaxsimB changes the random number seed
##D nsim<-50
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D Ndose<-length(doselev)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-4.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
##D 		sigmaup=30,edDF=5)
##D mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)
##D 
##D D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)
##D predict(D1[1],dose=c(75,125))
## End(Not run)



cleanEx()
nameEx("predict.emaxsimobj")
### * predict.emaxsimobj

flush(stderr()); flush(stdout())

### Name: predict.emaxsimobj
### Title: Mean response and SE for specified doses for a simulated
###   emaxsimobj object
### Aliases: predict.emaxsimobj
### Keywords: nonlinear

### ** Examples


## Not run: 
##D ## emaxsim changes the random number seed
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop.parm<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop.parm)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
##D D1 <- emaxsim(nsim,gen.parm)
##D d10<-D1[10]
##D predict(d10,c(75,150))
## End(Not run)
## Don't show: 
## emaxsim changes the random number seed
nsim<-3
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop.parm<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop.parm)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen.parm,nproc=1)
d3<-D1[3]
predict(d3,c(75,150))
## End(Don't show)



cleanEx()
nameEx("predict.fitEmax")
### * predict.fitEmax

flush(stderr()); flush(stdout())

### Name: predict.fitEmax
### Title: Estimated mean/proportion and confidence intervals derived from
###   the maximum likelihood fit of a 3- or 4- parameter Emax model.
### Aliases: predict.fitEmax
### Keywords: nonlinear

### ** Examples


## Not run: 
##D ## this example changes the random number seed
##D doselev<-c(0,5,25,50,100,350)
##D n<-c(78,81,81,81,77,80)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-8.0
##D pop.parm<-c(log(ed50),emax,e0)    
##D dose<-rep(doselev,n)
##D meanlev<-emaxfun(dose,pop.parm)  
##D 
##D y<-rnorm(sum(n),meanlev,sdy)
##D 
##D testout<-fitEmax(y,dose,modType=4)
##D predout<-predict(testout,dosevec=c(20,80),int=1)
## End(Not run)
## Don't show: 
## this example changes the random number seed
doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop.parm<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop.parm)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)
predout<-predict(testout,dosevec=c(20,80),int=1)
## End(Don't show)



cleanEx()
nameEx("predict.fitEmaxB")
### * predict.fitEmaxB

flush(stderr()); flush(stdout())

### Name: predict.fitEmaxB
### Title: Estimated mean and posterior intervals derived from a Bayesian
###   hyperbolic or sigmiodial Emax model.
### Aliases: predict.fitEmaxB

### ** Examples


## Not run: 
##D data("examples14")
##D exdat<-examples14[[6]]
##D 
##D prior<-prior.control(epmu=0,epsd=10,emaxmu=0,emaxsd=10,p50=0.25,
##D 				sigmalow=0.01,sigmaup=3)
##D mcmc<-mcmc.control(chains=3)
##D 
##D fitout<-fitEmaxB(exdat$y,exdat$dose,prior,modType=4,prot=exdat$prot,
##D 				count=exdat$nsize,msSat=(exdat$sd)^2,mcmc=mcmc)
##D 
##D predout<-predict(fitout,dosevec=sort(unique(exdat$dose)))
## End(Not run)



cleanEx()
nameEx("print.emaxsim")
### * print.emaxsim

flush(stderr()); flush(stdout())

### Name: print.emaxsim
### Title: Print simulation output from emaxsim
### Aliases: print.emaxsim
### Keywords: nonlinear

### ** Examples



## Not run: 
##D ## emaxsim changes the random number seed
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop.parm<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop.parm)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
##D D1 <- emaxsim(nsim,gen.parm)
##D 
##D print(D1,c(31,50),digits=2,id=4)
##D 
##D print(D1,c(1,20))
##D 
##D D1  ### implicitly calls print with default parameter settings
## End(Not run)
## Don't show: 
## emaxsim changes the random number seed
nsim<-3
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop.parm<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop.parm)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen.parm,nproc=1)

print(D1,c(2,3),digits=2,id=4)
## End(Don't show)



cleanEx()
nameEx("print.emaxsimB")
### * print.emaxsimB

flush(stderr()); flush(stdout())

### Name: print.emaxsimB
### Title: Print simulation output from emaxsimB
### Aliases: print.emaxsimB
### Keywords: nonlinear

### ** Examples


## Not run: 
##D ## emaxsimB changes the random number seed
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D Ndose<-length(doselev)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-4.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
##D 		sigmaup=30,edDF=5)
##D mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)
##D 
##D D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)
##D 
##D print(D1)
## End(Not run)



cleanEx()
nameEx("print.emaxsimobj")
### * print.emaxsimobj

flush(stderr()); flush(stdout())

### Name: print.emaxsimobj
### Title: Print a data set generated by emaxsim
### Aliases: print.emaxsimobj
### Keywords: nonlinear

### ** Examples


## Not run: 
##D 
##D save.seed<-.Random.seed
##D set.seed(12357)
##D 
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
##D D1 <- emaxsim(nsim,gen.parm)
##D e49<-D1[49]
##D 
##D e49
##D 
##D print(e49,c(101,200))  
##D 
##D .Random.seed<-save.seed
## End(Not run)



cleanEx()
nameEx("selEstan")
### * selEstan

flush(stderr()); flush(stdout())

### Name: selEstan
### Title: Select a pre-compiled 'rstan' Emax model
### Aliases: selEstan

### ** Examples

## Not run: 
##D 	estan<-selEstan()
## End(Not run)



cleanEx()
nameEx("showStanModels")
### * showStanModels

flush(stderr()); flush(stdout())

### Name: showStanModels
### Title: Display 'STAN' model code.
### Aliases: showStanModels

### ** Examples

## Not run: 
##D showStanModels()
## End(Not run)



cleanEx()
nameEx("sigma")
### * sigma

flush(stderr()); flush(stdout())

### Name: sigmaEmax
### Title: Extract Emax model residual SD estimates
### Aliases: sigma.fitEmax sigma.fitEmaxB sigma.emaxsim sigma.emaxsimB
### Keywords: nonlinear

### ** Examples

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)
sigma(testout)



cleanEx()
nameEx("startEmax")
### * startEmax

flush(stderr()); flush(stdout())

### Name: startEmax
### Title: Compute starting parameter values for the 3- or 4- Emax model.
### Aliases: startEmax
### Keywords: nonlinear

### ** Examples


data("examples14")
exdat<-examples14[[6]]
startEmax(exdat$y,exdat$dose)




cleanEx()
nameEx("summary.emaxsim")
### * summary.emaxsim

flush(stderr()); flush(stdout())

### Name: summary.emaxsim
### Title: Summary of output of emaxsim
### Aliases: summary.emaxsim
### Keywords: nonlinear

### ** Examples


## Not run: 
##D ## emaxsim changes the random number seed
##D nsim<-50
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop.parm<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop.parm)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
##D D1 <- emaxsim(nsim,gen.parm)
##D summary(D1,testalph=0.05,clev=0.95)
## End(Not run)
## Don't show: 
## emaxsim changes the random number seed
nsim<-3
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop.parm<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop.parm)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen.parm,nproc=1)
summary(D1,testalph=0.05,clev=0.95)
## End(Don't show)



cleanEx()
nameEx("summary.emaxsimB")
### * summary.emaxsimB

flush(stderr()); flush(stdout())

### Name: summary.emaxsimB
### Title: Summary of output of emaxsimB
### Aliases: summary.emaxsimB
### Keywords: nonlinear

### ** Examples


## Not run: 
##D 
##D ## emaxsimB changes the random number seed
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D Ndose<-length(doselev)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-4.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
##D 		sigmaup=30,edDF=5)
##D mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)
##D 
##D D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)
##D 
##D 
##D summary(D1,testalph=0.05,clev='0.95')
## End(Not run)



cleanEx()
nameEx("summary.emaxsimBobj")
### * summary.emaxsimBobj

flush(stderr()); flush(stdout())

### Name: summary.emaxsimBobj
### Title: Summarize Emax fit to a data set generated by emaxsimB
### Aliases: summary.emaxsimBobj
### Keywords: nonlinear

### ** Examples


## Not run: 
##D 
##D ## emaxsimB changes the random number seed
##D nsim<-50
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D Ndose<-length(doselev)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-4.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D prior<-prior.control(epmu=0,epsd=30,emaxmu=0,emaxsd=30,p50=50,sigmalow=0.1,
##D 		sigmaup=30,edDF=5)
##D mcmc<-mcmc.control(chains=1,warmup=500,iter=5000,seed=53453,propInit=0.15,adapt_delta = 0.95)
##D 
##D D1 <- emaxsimB(nsim,gen, prior, modType=3,mcmc=mcmc,check=FALSE)
##D summary(D1[1])
## End(Not run)



cleanEx()
nameEx("summary.emaxsimobj")
### * summary.emaxsimobj

flush(stderr()); flush(stdout())

### Name: summary.emaxsimobj
### Title: Summarize Emax fit to a data set generated by emaxsim
### Aliases: summary.emaxsimobj
### Keywords: nonlinear

### ** Examples


## emaxsim changes the random number seed
nsim<-3
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen.parm,nproc=1)
e3<-D1[3]

summary(e3)  



cleanEx()
nameEx("targetBeta")
### * targetBeta

flush(stderr()); flush(stdout())

### Name: targetBeta
### Title: Find a scaled Beta distribution matching specified probabilities
### Aliases: targetBeta
### Keywords: nonlinear

### ** Examples

### set quartiles at .15 and 1.0 for a beta distribution on (0,3)
targetBeta(minval=.15,pminV=0.25,pmaxV=0.75,maxval=1.0,upB=3)



cleanEx()
nameEx("targetCI")
### * targetCI

flush(stderr()); flush(stdout())

### Name: targetCI
### Title: Compute the dose with confidence interval exceeding a target
###   change from placebo for each simulated example in an emaxsim object.
### Aliases: targetCI
### Keywords: nonlinear

### ** Examples

	## Not run: 
##D 
##D 		# emaxsim changes the random number seed
##D 		nsim<-100
##D 		doselev<-c(0,5,25,50,100)
##D 		n<-c(78,81,81,81,77)
##D 
##D 		### population parameters for simulation
##D 		e0<-2.465375 
##D 		ed50<-67.481113 
##D 		emax<-15.127726
##D 		sdy<-7.967897
##D 		pop<-c(log(ed50),emax,e0)    
##D 		meanlev<-emaxfun(doselev,pop)  
##D 
##D 		###FixedMean is specialized constructor function for emaxsim
##D 		gen.parm<-FixedMean(n,doselev,meanlev,sdy)  
##D 
##D 		D1 <- emaxsim(nsim,gen.parm,modType=3)
##D 
##D 		target<-6
##D 		tD<- ( (target*ed50)/(emax-target) )
##D 		selectedDose<-targetCI(D1,target,dgrid=c(1:100)+0.5,cilev=0.80,high=TRUE)
##D 	
## End(Not run)
	## Don't show: 

		# emaxsim changes the random number seed
		nsim<-3
		doselev<-c(0,5,25,50,100)
		n<-c(78,81,81,81,77)

		### population parameters for simulation
		e0<-2.465375 
		ed50<-67.481113 
		emax<-15.127726
		sdy<-7.967897
		pop<-c(log(ed50),emax,e0)    
		meanlev<-emaxfun(doselev,pop)  

		###FixedMean is specialized constructor function for emaxsim
		gen.parm<-FixedMean(n,doselev,meanlev,sdy)  

		D1 <- emaxsim(nsim,gen.parm,modType=3,nproc=1)

		target<-6
		tD<- ( (target*ed50)/(emax-target) )
		selectedDose<-targetCI(D1,target,dgrid=c(1:100)+0.5,cilev=0.80,high=TRUE)
	
## End(Don't show)



cleanEx()
nameEx("targetD")
### * targetD

flush(stderr()); flush(stdout())

### Name: targetD
### Title: Compute the MLE (and its SE) of the dose achieving a specified
###   target improvement from placebo.
### Aliases: targetD
### Keywords: nonlinear

### ** Examples

    ## Not run: 
##D 
##D 		## emaxsim changes the random number seed
##D     doselev<-c(0,5,25,50,100,250)
##D     n<-c(78,81,81,81,77,80)
##D     dose<-rep(doselev,n)
##D 
##D     ### population parameters for simulation
##D     e0<-2.465375 
##D     ed50<-67.481113 
##D     emax<-15.127726
##D     sdy<-7.967897
##D     pop<-c(led50=log(ed50),emax=emax,e0=e0)    
##D     meanresp<-emaxfun(dose,pop)  
##D     y<-rnorm(sum(n),meanresp,sdy)
##D     nls.fit<-nls(y ~ e0 + (emax * dose)/(dose + exp(led50)), 
##D                               start = pop, control = nls.control(
##D                               maxiter = 100),trace=TRUE,na.action=na.omit)
##D 
##D     targetD(nls.fit,10,modType=3)
##D 
##D     ###
##D     ### apply targetD to an emaxsim object
##D     ###
##D     nsim<-50
##D     sdy<-25
##D     gen.parm<-FixedMean(n,doselev,emaxfun(doselev,pop),sdy)  
##D     D4 <- emaxsim(nsim,gen.parm,modType=4)
##D     summary(D4,testalph=0.05)
##D 
##D     out<-NULL
##D     for(i in 1:nsim){
##D         out<-rbind(out,targetD(D4[i],target=4))
##D     }
## End(Not run)
    ## Don't show: 

		## emaxsim changes the random number seed
    doselev<-c(0,5,25,50,100,250)
    n<-c(78,81,81,81,77,80)
    dose<-rep(doselev,n)

    ### population parameters for simulation
    e0<-2.465375 
    ed50<-67.481113 
    emax<-15.127726
    sdy<-7.967897
    pop<-c(led50=log(ed50),emax=emax,e0=e0)    
    meanresp<-emaxfun(dose,pop)  
    y<-rnorm(sum(n),meanresp,sdy)
    nls.fit<-nls(y ~ e0 + (emax * dose)/(dose + exp(led50)), 
                              start = pop, control = nls.control(
                              maxiter = 100),trace=TRUE,na.action=na.omit)

    targetD(nls.fit,10,modType=3)
## End(Don't show)



cleanEx()
nameEx("update.emaxsimobj")
### * update.emaxsimobj

flush(stderr()); flush(stdout())

### Name: update.emaxsimobj
### Title: Update estimation in a data set generated by emaxsim
### Aliases: update.emaxsimobj
### Keywords: nonlinear

### ** Examples


## Not run: 
##D 
##D ## emaxsim changes the random number seed
##D nsim<-50
##D idmax<-5
##D doselev<-c(0,5,25,50,100)
##D n<-c(78,81,81,81,77)
##D 
##D ### population parameters for simulation
##D e0<-2.465375 
##D ed50<-67.481113 
##D emax<-15.127726
##D sdy<-7.967897
##D pop<-c(log(ed50),emax,e0)    
##D meanlev<-emaxfun(doselev,pop)  
##D 
##D ###FixedMean is specialized constructor function for emaxsim
##D gen<-FixedMean(n,doselev,meanlev,sdy)  
##D D1 <- emaxsim(nsim,gen)
##D e49<-D1[49]
##D 
##D #### re-try estimation starting at the population value
##D e49u<- update(e49,pop)
## End(Not run)
## Don't show: 
## emaxsim changes the random number seed
nsim<-3
doselev<-c(0,5,25,50,100)
n<-c(78,81,81,81,77)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-7.967897
pop<-c(log(ed50),emax,e0)    
meanlev<-emaxfun(doselev,pop)  

###FixedMean is specialized constructor function for emaxsim
gen<-FixedMean(n,doselev,meanlev,sdy)  
D1 <- emaxsim(nsim,gen,nproc=1)
e3<-D1[3]

#### re-try estimation starting at the population value
e3u<- update(e3,pop)
## End(Don't show)



cleanEx()
nameEx("vcov")
### * vcov

flush(stderr()); flush(stdout())

### Name: vcovEmax
### Title: Extract Emax model variance-covariance matrix for ML estimates
### Aliases: vcov.fitEmax vcov.emaxsim
### Keywords: nonlinear

### ** Examples

doselev<-c(0,5,25,50,100,350)
n<-c(78,81,81,81,77,80)

### population parameters for simulation
e0<-2.465375 
ed50<-67.481113 
emax<-15.127726
sdy<-8.0
pop<-c(log(ed50),emax,e0)    
dose<-rep(doselev,n)
meanlev<-emaxfun(dose,pop)  

y<-rnorm(sum(n),meanlev,sdy)

testout<-fitEmax(y,dose,modType=4)
vcov(testout)



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
