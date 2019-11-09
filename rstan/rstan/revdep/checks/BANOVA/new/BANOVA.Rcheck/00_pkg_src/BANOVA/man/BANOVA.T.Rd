\name{BANOVA.T}
\alias{BANOVA.T}
\alias{predict.BANOVA.T}
\alias{print.BANOVA.T}
\alias{summary.BANOVA.T}
\title{Estimation of BANOVA with T-distributin of the dependent variable}
\description{
\code{BANOVA.T} implements a Hierarchical Bayesian ANOVA for linear models with T-distributed response. 
}
\usage{
BANOVA.T(l1_formula = "NA", l2_formula = "NA", data, id, l1_hyper = c(1, 1, 1), 
l2_hyper = c(1, 1, 0.0001), burnin = 5000, sample = 2000, thin = 10, 
adapt = 0, conv_speedup = F, jags = runjags.getOption('jagspath'))
\method{summary}{BANOVA.T}(object, ...)
\method{predict}{BANOVA.T}(object, newdata = NULL,...)
\method{print}{BANOVA.T}(x, ...)
}
\arguments{
     \item{l1_formula}{formula for level 1 e.g. 'Y~X1+X2'}

  \item{l2_formula}{formula for level 2 e.g. '~Z1+Z2',
  response variable must not be included, if missing, the single level model will be generated}

  \item{data}{a data.frame in long format including all features in level 1 and level 2(covariates and categorical factors) and responses}

  \item{id}{subject ID of each response unit}
  
  \item{l1_hyper}{level 1 hyperparameters, c(\eqn{\alpha, \beta, \lambda}) for two-level models and c(\eqn{\alpha, \beta, \lambda, \sigma_p}) for single level models, default c(1,1,1)}
  
  \item{l2_hyper}{level 2 hyperparameters, c(a, b, \eqn{\gamma}), default c(1,1,0.0001)}

  \item{burnin}{the number of burn in draws in the MCMC algorithm, default 5000}

  \item{sample}{target samples in the MCMC algorithm after thinning, default 2000}

  \item{thin}{the number of samples in the MCMC algorithm that needs to be thinned, default 10}
  
  \item{adapt}{the number of adaptive iterations, default 0 (see \link[runjags]{run.jags})}
  
  \item{conv_speedup}{whether to speedup convergence, default F}
  
  \item{jags}{the system call or path for activating 'JAGS'. Default calls findjags() to attempt to locate 'JAGS' on your system}
  
  \item{object}{object of class \code{BANOVA.T} (returned by \code{BANOVA.T})}

  \item{newdata}{test data, either a matrix, vector or a
  data frame. It must have the same format with the original data (the same column number)}
  
  \item{x}{object of class \code{BANOVA.T} (returned by \code{BANOVA.T})}
  
  \item{\dots}{additional arguments,currently ignored}
}

\details{
Level 1 model: \cr
  \eqn{y_i} {~} \eqn{t(\nu, \eta_i,\sigma^{-2})} \cr
  where \eqn{\eta_i = \sum_{p = 0}^{P}\sum_{j=1}^{J_p}X_{i,j}^p\beta_{j,s_i}^p}, \eqn{s_i} is the subject id of response \eqn{i}, see \code{\link{BANOVA-package}}. The hyper parameters: \eqn{\nu} is the degree of freedom, \eqn{\nu} {~} {Piosson(\eqn{\lambda})} and \eqn{\sigma} is the scale parameter, \eqn{\sigma^{-2}} {~} {Gamma(\eqn{\alpha, \beta})}.
  
}
\value{
  \code{BANOVA.T} returns an object of class \code{"BANOVA.T"}. The returned object is a list containing:
  
  \item{anova.table}{table of effect sizes \code{\link{BAnova}}}
  \item{coef.tables}{table of estimated coefficients}
  \item{pvalue.table}{table of p-values \code{\link{table.pvalues}}}
  \item{dMatrice}{design matrices at level 1 and level 2}
  \item{samples_l2_param}{posterior samples of level 2 parameters}
  \item{data}{original data.frame}
  \item{mf1}{model.frame of level 1}
  \item{mf2}{model.frame of level 2}
  \item{JAGSmodel}{'JAGS' model}
}
\examples{
\donttest{
# Use the ipadstudy data set
data(ipadstudy)
res <- BANOVA.T(attitude~1, ~owner + age + gender + selfbrand*conspic, ipadstudy,
ipadstudy$id, burnin = 5000, sample = 2000, thin = 10)
summary(res)

# or use BANOVA.run based on 'Stan'
res19 <- BANOVA.run(attitude~1, ~owner + age + gender + selfbrand*conspic, 
data = ipadstudy, model_name = 'T', id = ipadstudy$id, iter = 100, 
thin = 1, chains = 2)
}
}

