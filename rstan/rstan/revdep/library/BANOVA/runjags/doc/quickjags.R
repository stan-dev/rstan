## ---- echo = FALSE, include=FALSE----------------------------------------
library(runjags)
if(require(rjags)){
	load.module('dic')
}
library(parallel)

## ------------------------------------------------------------------------
testjags()

## ---- echo = FALSE, include=FALSE----------------------------------------
runjags.options(inits.warning=FALSE, rng.warning=FALSE, blockignore.warning=FALSE, blockcombine.warning=FALSE, silent.jags=TRUE, silent.runjags=TRUE)

## ------------------------------------------------------------------------
filestring <- "
Poisson model...

model {
for (i in 1:6) {
for (j in 1:3) {
y[i,j] ~ dpois(mu[i])
}
log(mu[i]) <- alpha + beta*log(x[i] + 10) + gamma*x[i]
}
for (i in 1:6) {
y.pred[i] ~ dpois(mu[i])
}
alpha ~ dnorm(0, 0.0001)
beta ~ dnorm(0, 0.0001)
gamma ~ dnorm(0, 0.0001)
}

Data{
list(y = structure(.Data = c(15,21,29,16,18,21,16,26,33,
27,41,60,33,38,41,20,27,42),
.Dim = c(6, 3)),
x = c(0, 10, 33, 100, 333, 1000))
}

Inits{
list(alpha = 0, beta = 1, gamma = -0.1)
}

Inits{
list(alpha = 10, beta = 0, gamma = 0.1)
}
"

## ---- results='hide'-----------------------------------------------------
results <- run.jags(filestring, monitor=c('alpha','beta','gamma'))

## ------------------------------------------------------------------------
results

## ---- results='hide', include=FALSE--------------------------------------
plot(results, layout=c(3,4), file='plots.pdf', height=10, width=8)

## ---- results='hide'-----------------------------------------------------
plot(results, layout=c(3,4))

## ---- results='hide'-----------------------------------------------------
resultswithglm <- run.jags(filestring, monitor=c('alpha','beta','gamma'), modules='glm')

## ------------------------------------------------------------------------
resultswithglm

## ---- results='hide'-----------------------------------------------------
# Simulate the data:
set.seed(1)
X <- 1:100
Y <- rnorm(length(X), 2*X + 10, 1)
N <- length(X)

model <- "
	model {
        for(i in 1 : N){ #data# N
        Y[i] ~ dnorm(true.y[i], precision) #data# Y
        true.y[i] <- (coef * X[i]) + int #data# X
	}
  	coef ~ dunif(-1000,1000)
  	int ~ dunif(-1000,1000)
  	precision ~ dexp(1)
  	#inits# coef, int, precision
  	#monitor# coef, int, precision
}"
# A function to return initial values for each chain:
coef <- function(chain) return(switch(chain, "1"= -10, "2"= 10))
int <- function(chain) return(switch(chain, "1"= -10, "2"= 10))
precision <- function(chain) return(switch(chain, "1"= 0.01, "2"= 100))

# Run the simulation:
results <- run.jags(model, n.chains = 2)

## ---- results='hide'-----------------------------------------------------
results <- extend.jags(results, sample=5000)

## ---- results='hide'-----------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
D9 <- data.frame(weight, group)

# The JAGS equivalent:
model <- template.jags(weight ~ group, D9, n.chains=2, family='gaussian')

## ---- results='hide'-----------------------------------------------------
JAGS.D9 <- run.jags(model)
lm.D9 <- lm(weight ~ group, data=D9)

## ------------------------------------------------------------------------
JAGS.D9
summary(lm.D9)

## ---- results='hide'-----------------------------------------------------
JAGS.D9 <- run.jags(model, mutate=list(prec2sd, 'precision'))

## ------------------------------------------------------------------------
summary(JAGS.D9, vars=c('regression_precision.sd', 'intercept', 'group_effect'))

## ------------------------------------------------------------------------
summary(residuals(lm.D9) - residuals(JAGS.D9, output='mean'))

## ------------------------------------------------------------------------
extract(JAGS.D9, what='samplers')

## ---- eval=FALSE---------------------------------------------------------
#  ?runjags.options

## ------------------------------------------------------------------------
citation('runjags')

## ------------------------------------------------------------------------
sessionInfo()

