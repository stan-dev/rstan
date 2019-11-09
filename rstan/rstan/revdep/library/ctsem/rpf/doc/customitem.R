## ---- message=FALSE------------------------------------------------------
library(rpf)

lmp.item<-rpf.lmp(k=2) # create item w/ 5th order polynomial
par<-c(.69,.71,-.5,-8.48,.52,-3.32) # item parameters
theta<-seq(-3,3,.1) # grid for latent trait

## Test the traceline or "prob" C++ function
P<-rpf.prob(lmp.item, par, theta)

## Prettier plots than this are of course possible
plot(theta, P[2,], type="l", ylim=c(0,1), xlab="Theta", ylab="P(Theta)")


## ------------------------------------------------------------------------
## Derivatives of negative log-likelihood at arbitrary point with arbitrary weights
## Rounding only for easy reading for tutorial
round(rpf.dLL(lmp.item, par, theta[1], weight=c(5,7)),2)

## ------------------------------------------------------------------------
## Analytical derivatives "deriv1" followed by "deriv2"
rpf.dTheta(lmp.item, par, where=-.5, dir=1)

## Numerical derivatives
library(numDeriv)
dTheta.wrap<-function(theta, spec, par, cat=1){
  rpf.prob(spec, par, theta)[cat]
}

## should match first element of gradient from rpf.dTheta
grad(dTheta.wrap, -.5, spec=lmp.item, par=par)

## should match first element of hessian from rpf.dTheta
hessian(dTheta.wrap, -.5, spec=lmp.item, par=par)


