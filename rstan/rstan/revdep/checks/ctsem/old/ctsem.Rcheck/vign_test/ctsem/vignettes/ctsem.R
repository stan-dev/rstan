## ----setup, include = FALSE, cache = FALSE, echo = FALSE----------------------
library('ctsem')
library('knitr')
render_sweave()
set.seed(22)
knit_hooks$set(crop = hook_pdfcrop)
opts_chunk$set(fig.path = 'figures/plots-', warning = FALSE, fig.align = 'center', width.cutoff = 80, fig.show = 'hold', eval = TRUE, echo = TRUE, message = FALSE, background = "white", prompt = FALSE, highlight = FALSE, comment = NA, tidy = FALSE, out.truncate = 80)
options(replace.assign = TRUE, width = 80, prompt = "R> ", scipen = 12, digits = 3,crop=TRUE)


# setwd('C:\\Users\\driver\\Dropbox\\MPIB\\CT-SEM\\manual') #set this working directory!
Sys.setenv(TEXINPUTS = getwd(),
  BIBINPUTS = getwd(),
  BSTINPUTS = getwd())

## ----eval=FALSE---------------------------------------------------------------
#  data('datastructure')
#  datastructure
#  semModel<-ctModel(n.latent=2, n.manifest=3, TRAITVAR='auto',
#    n.TIpred=2, n.TDpred=1, Tpoints=3,
#    LAMBDA=matrix(c(1,'lambda21', 0, 0,1,0),nrow=3))
#  semFit<-ctFit(datastructure, semModel, nofit=TRUE)
#  semFit$mxobj$A$labels
#  semFit$mxobj$S$labels
#  semFit$mxobj$M$labels
#  semFit$mxobj$F$values

## ----install, echo = TRUE, eval = FALSE---------------------------------------
#  install.packages("ctsem")
#  library("ctsem")

## ----wideformat, echo = FALSE, out.truncate=100, width.cutoff=100-------------
options(width = 100)
data('datastructure')
datastructure
options(width = 80)

## ----longformat, include = TRUE, cache = FALSE, echo = FALSE, results = 'markup'----
data('longexample')
head(longexample, 7)

## ----longformatconversion, include = TRUE, cache = FALSE, echo = TRUE, results = 'markup'----
data("longexample")
wideexample <- ctLongToWide(datalong = longexample, id = "id", 
  time = "time", manifestNames = c("Y1", "Y2", "Y3"), 
  TDpredNames = "TD1", TIpredNames = c("TI1", "TI2"))
wide <- ctIntervalise(datawide = wideexample, Tpoints = 4, n.manifest = 3, 
  n.TDpred = 1, n.TIpred = 2, manifestNames = c("Y1", "Y2", "Y3"), 
  TDpredNames = "TD1", TIpredNames = c("TI1", "TI2") )

## ----simplemodel, include = TRUE, echo = TRUE, results = 'hide'---------------
examplemodel <- ctModel(n.latent = 2, n.manifest = 2, Tpoints = 3, 
  LAMBDA = diag(2))

## ----example1ctfit, include = TRUE, cache = TRUE, echo = TRUE, results = 'hide'----
data("ctExample1")
example1model <- ctModel(n.latent = 2, n.manifest = 2, Tpoints = 6, 
  manifestNames = c("LeisureTime", "Happiness"), 
  latentNames = c("LeisureTime", "Happiness"), LAMBDA = diag(2))
example1fit <- ctFit(dat = ctExample1, ctmodelobj = example1model)

## ----example1ctfittable, include = TRUE, echo = TRUE--------------------------
summary(example1fit, verbose = TRUE)["discreteDRIFTstd"]

## ----expm,eval=FALSE----------------------------------------------------------
#  expm(summary(example1fit)$DRIFT * 2.5)

## ----example1testing, cache = TRUE, echo = TRUE-------------------------------
testmodel <- example1model
testmodel$DRIFT[1, 2] <- 0
testfit <- ctFit(dat = ctExample1, ctmodelobj = testmodel)

## ----mxcompare----------------------------------------------------------------
mxCompare(example1fit$mxobj, testfit$mxobj)

## ----confidenceintervals, cache = TRUE, echo = 1------------------------------
example1cifit <- ctCI(example1fit, confidenceintervals = "DRIFT")
summary(example1cifit)$omxsummary$CI

## ----example2fit, cache = TRUE------------------------------------------------
data("ctExample1")
traitmodel <- ctModel(n.manifest = 2, n.latent = 2, Tpoints = 6, 
  LAMBDA = diag(2), manifestNames = c("LeisureTime", "Happiness"), 
  latentNames = c("LeisureTime", "Happiness"), TRAITVAR = "auto")
traitfit <- ctFit(dat = ctExample1, ctmodelobj = traitmodel)

## ----traitparamplot, include = TRUE, cache = FALSE, echo = FALSE, results = 'hide', fig.height = 4----
par(mfrow = c(2, 2))
par(mar = c(3, 4, 2, 2),mgp=c(1.5,.5,0),lwd=2)
plot(example1fit, wait = FALSE, max.time = 20, mean = FALSE, withinVariance = FALSE,randomImpulse = FALSE,
  experimentalImpulse = FALSE)
plot(traitfit, wait = FALSE, max.time = 20, mean = FALSE, withinVariance = FALSE,randomImpulse = FALSE,
  experimentalImpulse = FALSE)

## ----example1TIpred, include = TRUE, cache = TRUE, echo = TRUE, results = 'hide'----
data("ctExample1TIpred")
tipredmodel <- ctModel(n.manifest = 2, n.latent = 2, n.TIpred = 1,
  manifestNames = c("LeisureTime", "Happiness"),
  latentNames = c("LeisureTime", "Happiness"),
  TIpredNames = "NumFriends",
 Tpoints = 6, LAMBDA = diag(2), TRAITVAR = "auto")
tipredfit <- ctFit(dat = ctExample1TIpred, ctmodelobj = tipredmodel)

summary(tipredfit, verbose = TRUE)["TIPREDEFFECT"]
summary(tipredfit, verbose = TRUE)["discreteTIPREDEFFECT"]
summary(tipredfit, verbose = TRUE)["asymTIPREDEFFECT"]
summary(tipredfit, verbose = TRUE)["addedTIPREDVAR"]

## ----example1TIpredestimates, echo = FALSE, out.width = '4cm', out.height = '4cm'----
summary(tipredfit, verbose = TRUE)["TIPREDEFFECT"]
cat("\n")
summary(tipredfit, verbose = TRUE)["discreteTIPREDEFFECT"]

## ----example1TIpredestimates2, echo = FALSE, out.width = '4cm', out.height = '4cm'----
summary(tipredfit, verbose = TRUE)["asymTIPREDEFFECT"]
cat("\n")
summary(tipredfit, verbose = TRUE)["addedTIPREDVAR"]

## ----tdpreddemo, echo = FALSE, eval = TRUE, fig.height = 3--------------------
par(mar=c(3.0,2.5,1.5,.5)+.1,mgp=c(1.5,.5,0),lwd=2)
Tpoints = 19
testm <- ctModel(Tpoints = Tpoints, n.latent = 1, 
  n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.5, 1), TDPREDEFFECT = diag(1, 1), 
  CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(100, 1), 
  TDPREDMEANS = matrix(c(0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ncol = 1, nrow = (Tpoints)))
testd <- ctGenerate(testm, n.subjects = 100, burnin = 30)

testm <- ctModel(Tpoints = Tpoints, n.latent = 1, n.TDpred = 1, 
  n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.001, 1), TDPREDEFFECT = diag(1, 1), 
  CINT = diag(4, 1), 
  T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(1, 1), 
  TDPREDMEANS = matrix(c(0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ncol = 1, nrow = (Tpoints)))
testdpure <- ctGenerate(testm, n.subjects = 1, burnin = 30)

par(mfrow = c(1, 2), cex = .7)
plot(0:(Tpoints-1), testdpure[, 1:Tpoints], type = 'l', ylim = c(min(testd[, 1:Tpoints]), max(testd[, 1:Tpoints])),
  ylab = 'Dependent variable', xlab = 'Time', lwd = 3, main = 'Impulse predictor')
for(i in 1:5){
  points(0:(Tpoints-1), testd[i, 1:Tpoints], col = 1+i, type = 'b')
}

Tpoints = 19
testm <- ctModel(Tpoints = Tpoints, n.latent = 1, 
  n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.5, 1), TDPREDEFFECT = diag(1.6, 1), 
  CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(100, 1), 
  TDPREDMEANS = matrix(c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), ncol = 1, nrow = (19)))
testd <- ctGenerate(testm, n.subjects = 100, burnin = 30)

testm <- ctModel(Tpoints = Tpoints, n.latent = 1, 
  n.TDpred = 1, n.manifest = 1, LAMBDA = diag(1), DRIFT = diag(-.3, 1),
  DIFFUSION = diag(.001, 1), TDPREDEFFECT = diag(1.6, 1), 
  CINT = diag(4, 1), T0MEANS = matrix(0, ncol = 1, nrow = 1),
  T0VAR = diag(1, 1), 
  TDPREDMEANS = matrix(c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), ncol = 1, nrow = (19)))
testdpure <- ctGenerate(testm, n.subjects = 1, burnin = 30)

plot(0:(Tpoints-1), testdpure[, 1:Tpoints], type = 'l', ylim = c(min(testd[, 1:Tpoints]), max(testd[, 1:Tpoints])), ylab = 'Dependent variable', xlab = 'Time', lwd = 3, main = 'Level predictor')
for(i in 1:5){
  points(0:(Tpoints-1), testd[i, 1:Tpoints], col = 1+i, type = 'b')
}

## ----example2TDpred, include = TRUE, cache = TRUE, echo = TRUE, results = 'show', fig.height = 4, fig.width = 4, fig.align = 'center'----
data("ctExample2")
tdpredmodel <- ctModel(n.manifest = 2, n.latent = 2, n.TDpred = 1, 
  Tpoints = 8, manifestNames = c("LeisureTime", "Happiness"), 
  TDpredNames = "MoneyInt", latentNames = c("LeisureTime", "Happiness"),
  LAMBDA = diag(2), TRAITVAR = "auto")
tdpredfit <- ctFit(dat = ctExample2, ctmodelobj = tdpredmodel,
  stationary=c('T0VAR','T0TRAITEFFECT'))
summary(tdpredfit, verbose = TRUE)["TDPREDEFFECT"]


## ----example2TDpredlevel, include = TRUE, eval = TRUE, echo = TRUE, cache = TRUE----
data("ctExample2")
tdpredlevelmodel <- ctModel(n.manifest = 2, n.latent = 3, 
  n.TDpred = 1, 
  Tpoints = 8, manifestNames = c("LeisureTime", "Happiness"), 
  TDpredNames = "MoneyInt", 
  latentNames = c("LeisureTime", "Happiness", "MoneyIntLatent"),
  LAMBDA = matrix(c(1,0, 0,1, 0,0), ncol = 3), TRAITVAR = "auto")

tdpredlevelmodel$TRAITVAR[3, ] <- 0
tdpredlevelmodel$TRAITVAR[, 3] <- 0
tdpredlevelmodel$DIFFUSION[, 3] <- 0
tdpredlevelmodel$DIFFUSION[3, ] <- 0
tdpredlevelmodel$T0VAR[3, ] <- 0
tdpredlevelmodel$T0VAR[, 3] <- 0
tdpredlevelmodel$CINT[3] <- 0
tdpredlevelmodel$T0MEANS[3] <- 0
tdpredlevelmodel$TDPREDEFFECT[1:3, ] <- c(0,0,1)
tdpredlevelmodel$DRIFT[3, ] <- c(0,0,-.000001)

tdpredlevelfit <- ctFit(dat = ctExample2, 
  ctmodelobj = tdpredlevelmodel, 
  stationary=c('T0VAR','T0TRAITEFFECT'))

summary(tdpredlevelfit, verbose = TRUE)[c("DRIFT","TDPREDEFFECT")]

## ----timeseries, cache = TRUE, echo = TRUE------------------------------------
data("ctExample3")
model <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 100, 
  LAMBDA = matrix(c(1, "lambda2", "lambda3"), nrow = 3, ncol = 1))
fit <- ctFit(dat = ctExample3, ctmodelobj = model, objective = "Kalman",
  stationary = c("T0VAR"))

## ----multigroup, cache = TRUE, echo = TRUE------------------------------------
data("ctExample4")

basemodel <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 20,
  LAMBDA = matrix(c(1, "lambda2", "lambda3"), nrow = 3, ncol = 1))

freemodel <- basemodel
freemodel$LAMBDA[3, 1] <- "groupfree"
groups <- paste0("g", rep(1:2, each = 10))

multif <- ctMultigroupFit(dat = ctExample4, groupings = groups,
  ctmodelobj = basemodel, freemodel = freemodel)

## ----multigroupOutput, echo = FALSE-------------------------------------------
multif$mxobj$output$estimate[grep("lambda3", names(multif$mxobj$output$estimate))]

## ----dynamicresidualmodel, echo = TRUE, results='hide',fig.keep='none'--------
genm <- ctModel(Tpoints = 200, n.latent = 2, n.manifest = 1, 
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  DIFFUSION = matrix(c(0, 0, 0, 1), 2),
  MANIFESTVAR = t(chol(diag(.6,1))),
  DRIFT = matrix(c(0, -.1, 1, -.2), nrow = 2),   
  CINT = matrix(c(1, 0), nrow = 2))

data <- ctGenerate(genm, n.subjects = 1, burnin = 200)

ctIndplot(data, n.subjects = 1 , n.manifest = 1, Tpoints = 200)

model <- ctModel(Tpoints = 200, n.latent = 2, n.manifest = 1, 
  LAMBDA = matrix(c(1, 0), nrow = 1, ncol = 2),
  DIFFUSION = matrix(c(0, 0, 0, "diffusion"), 2),
  DRIFT = matrix(c(0, "regulation", 1, "diffusionAR"), nrow = 2))

fit <- ctFit(data, model, stationary = c("T0VAR"))

