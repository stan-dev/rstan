#' Diagnostics for ctsem importance sampling
#'
#' @param fit Output from ctStanFit when optimize=TRUE and isloops > 0
#'
#' @return Nothing. Plots convergence of parameter mean estimates from initial Hessian based distribution to final sampling distribution.
#' @export
#'
#' @examples
#' \dontrun{
#' #get data
#' sunspots<-sunspot.year
#' sunspots<-sunspots[50: (length(sunspots) - (1988-1924))]
#' id <- 1
#' time <- 1749:1924
#' datalong <- cbind(id, time, sunspots)
#'
#' #setup model
#' model <- ctModel(type='stanct', n.latent=2, n.manifest=1, 
#'  manifestNames='sunspots', 
#'  latentNames=c('ss_level', 'ss_velocity'),
#'   LAMBDA=matrix(c( -1, 'ma1' ), nrow=1, ncol=2),
#'   DRIFT=matrix(c(-.0001, 'a21', 1, 'a22'), nrow=2, ncol=2),
#'   MANIFESTMEANS=matrix(c('m1'), nrow=1, ncol=1),
#'   CINT=matrix(c(0, 0), nrow=2, ncol=1),
#'   T0VAR=matrix(c(1,0,0,1), nrow=2, ncol=2), #Because single subject
#'   DIFFUSION=matrix(c(0.0001, 0, 0, "diffusion"), ncol=2, nrow=2))
#' 
#' model$pars$indvarying<-FALSE #Because single subject
#' model$pars$transform[14]<- '(param)*5+44 ' #Because not mean centered
#' model$pars$transform[4]<-'-log(exp(-param*1.5)+1)' #To avoid multi modality
#'
#' #fit and plot importance sampling diagnostic
#' fit <- ctStanFit(datalong, model, chains=1,
#'   optimcontrol=list(isloops=5,finishsamples=500),optimize=TRUE)
#' isdiag(fit)
#' }
 
isdiag <- function(fit){
  iter=length(fit$stanfit$isdiags$cov)
  mcov <- fit$stanfit$isdiags$cov
  samplecov <- cov(fit$stanfit$rawposterior)
  means <- simplify2array(fit$stanfit$isdiags$means)
  means <- (means - means[,ncol(means)]) 
  means <- t(means / sqrt(diag(samplecov)))
  
  # smeans <- matrix(apply(means,2,function(x) x))),byrow=TRUE,ncol=iter)
matplot(means,type='l',main='Mean convergence',xlab='Sampling loop',ylab=' Z divergence relative to finish',xlim=c(0,iter*1.2))

legend('topright',bty='n',legend = paste0('par',1:ncol(means)),lty = 1:5,col=1:6,text.col=1:6,cex = .7)

sds <- simplify2array(mcov)
sds <- apply(sds,3,function(x) sqrt(diag(x)))
sds <- t((sds - sds[,iter]) / sqrt(diag(samplecov)))
matplot(sds,type='l',main='SD convergence',xlab='Sampling loop',ylab='Z divergence relative to finish',xlim=c(0,iter*1.2))

legend('topright',bty='n',legend = paste0('par',1:ncol(means)),lty = 1:5,col=1:6,text.col=1:6,cex = .7)




}
