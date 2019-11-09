require(ctsem)
require(testthat)

context("differenttraits")

test_that("time calc", {
set.seed(4)
Tpoints<-30
n.manifest=2
nsubjects=1000
n.latent=2

DRIFT=matrix(c(-.3, .2, 0, -0.5), byrow=TRUE, nrow=n.latent, ncol=n.latent)

genm=ctModel(Tpoints=Tpoints,
  n.latent=n.latent, n.manifest=n.manifest,
  LAMBDA=matrix(c(1, 0,0,1), nrow=n.manifest, ncol=n.latent),
  DRIFT=DRIFT,
  DIFFUSION=matrix(c(2, 0, 0, 1), byrow=TRUE, nrow=n.latent, ncol=n.latent),
  MANIFESTVAR=matrix(c(1, 0,0,.5), nrow=n.manifest, ncol=n.manifest),
  TRAITVAR=matrix(c(1,.5,0,.8),n.latent,n.latent))

cd=ctGenerate(ctmodelobj=genm, n.subjects=nsubjects, burnin=51, dtmean=1, 
  logdtsd=0,wide=TRUE)

wide=cd

long=ctWideToLong(datawide = cd,Tpoints = Tpoints,n.manifest = n.manifest)
long=ctDeintervalise(datalong = long)
long=long[-seq(3,length(long),3),]
wide=ctLongToWide(datalong = long,id='id',time='time',manifestNames= genm$manifestNames)

Tpoints=20
wide=ctIntervalise(datawide = wide,Tpoints = Tpoints,n.manifest = n.manifest)

mltrait<-ctModel(Tpoints=Tpoints,n.latent=n.latent,n.manifest=n.manifest,
  LAMBDA=diag(1,n.manifest),
  TRAITVAR='auto')

mmtrait<-ctModel(Tpoints=Tpoints,n.latent=n.latent,n.manifest=n.manifest,
  LAMBDA=diag(1,n.manifest),
  MANIFESTTRAITVAR='auto')

mptrait<-ctModel(Tpoints=Tpoints,n.latent=4,n.manifest=n.manifest,
  LAMBDA=matrix(c(1,0, 0,1, 0,0, 0,0),2,4),
  DRIFT=matrix(c(
    'dr11','dr12',1,0,
    'dr21','dr22',0,1,
    0,0,.0001,0,
    0,0,0,.0001),byrow=TRUE,4,4),
  DIFFUSION=matrix(c(
    'df11',0,0,0,
    'df21','df22',0,0,
    0,0,.0001,0,
    0,0,0,.0001),byrow=TRUE,4,4),
  T0MEANS=matrix(c('t1','t2',0,0),ncol=1))

fmlstrait=ctFit(dat = wide,ctmodelobj = mltrait,retryattempts = 5,stationary='T0TRAITEFFECT')
fmltrait=ctRefineTo(dat = wide,ctmodelobj = mltrait,retryattempts = 5,stationary='')

dfmlstrait=ctFit(dat= wide,ctmodelobj = mltrait,retryattempts = 5,discreteTime=TRUE,stationary='T0TRAITEFFECT')
dfmltrait=ctFit(dat = wide,ctmodelobj = mltrait,retryattempts = 5,discreteTime=TRUE,stationary='')

fmmtrait=ctFit(dat = wide,ctmodelobj = mmtrait,retryattempts = 5)
fmptrait=ctFit(dat = wide,ctmodelobj = mptrait,retryattempts = 5,stationary='')

# summary(fmlstrait)
# summary(fmmtrait)
# summary(fmptrait)
# summary(fmltrait)

#check traits using different fit approaches
expect_equal(rep(0,4),c(fmlstrait$mxobj$DRIFT$values-fmmtrait$mxobj$DRIFT$values),tolerance=1e-2)
expect_equal(rep(0,4),c(fmltrait$mxobj$DRIFT$values-fmptrait$mxobj$DRIFT$values[1:2,1:2]),tolerance=1e-2)

# expect_equal(rep(0,4),c(expm(fmlstrait$mxobj$DRIFT$values)-dfmltrait$mxobj$DRIFT$values),tolerance=1e-2)

#check DRIFT is reasonably estimated
expect_equal(rep(0,4),c(fmltrait$mxobj$DRIFT$values-DRIFT),tolerance=.1)

})
