\name{BANOVA.ordMultinomial}
\alias{BANOVA.ordMultinomial}
\alias{predict.BANOVA.ordMultinomial}
\alias{print.BANOVA.ordMultinomial}
\alias{summary.BANOVA.ordMultinomial}
\title{Estimation of BANOVA with a ordered Multinomial response variable}
\description{
\code{BANOVA.ordMultinomial} implements a Hierarchical Bayesian ANOVA for ordered multinomial responses, with
 a normal heterogeneity distribution.
}

\usage{
BANOVA.ordMultinomial(l1_formula = "NA",
  l2_formula = "NA", data, id, l1_hyper = c(0.0001, 100), 
  l2_hyper = c(1, 1, 0.0001, 100), burnin = 5000, 
  sample = 2000, thin = 10, adapt = 0, conv_speedup = F, 
  jags = runjags.getOption('jagspath'))
\method{summary}{BANOVA.ordMultinomial}(object, ...)
\method{predict}{BANOVA.ordMultinomial}(object, newdata = NULL,...)
\method{print}{BANOVA.ordMultinomial}(x, ...)
}
\arguments{
  \item{l1_formula}{formula for level 1 e.g. 'Y~X1+X2'}

  \item{l2_formula}{formula for level 2 e.g. '~Z1+Z2',
  response variable must not be included, if missing, the single level model will be generated}

  \item{data}{a data frame}

  \item{id}{subject ID of each response unit}
  
  \item{l1_hyper}{level 1 hyperparameters for single level models, default c(0.0001,100)}
  
  \item{l2_hyper}{level 2 hyperparameters, c(a, b, \eqn{\gamma}, d), default c(1,1,0.0001,100)}

  \item{burnin}{the number of burn in draws in the MCMC algorithm, default 5000}

  \item{sample}{target samples in the MCMC algorithm after thinning, default 2000}

  \item{thin}{the number of samples in the MCMC algorithm that needs to be thinned, default 10}
  
  \item{adapt}{the number of adaptive iterations, default 0 (see \link[runjags]{run.jags})}
  
  \item{conv_speedup}{whether to speedup convergence, default F}
  
  \item{jags}{the system call or path for activating 'JAGS'. Default calls findjags() to attempt to locate 'JAGS' on your system}
  
  \item{object}{object of class \code{BANOVA.ordMultinomial} (returned by \code{BANOVA.ordMultinomial})}

  \item{newdata}{test data, either a matrix, vector or a
  data frame. It must have the same format with the original data (the same column number)}
  
  \item{x}{object of class \code{BANOVA.ordMultinomial} (returned by \code{BANOVA.ordMultinomial})}
  
  \item{\dots}{additional arguments,currently ignored}
}

\details{
Level 1 model: \cr
  \eqn{y_i = 1}, if \eqn{l_i < 0} \cr
  \eqn{y_i = 2}, if \eqn{0 < l_i < c_2} \cr
  {...} \cr
  \eqn{y_i = \ell}, if \eqn{c_{\ell-1} < l_i < \infty} \cr
  \eqn{l_i = \eta_i + \epsilon_i}
  where \eqn{\epsilon_i} {~} logistic \eqn{(0,1)}, \eqn{c_\ell, (\ell = 2,...L-1)} are cut points, \eqn{c_\ell} {~} \eqn{N(0, \bar{\sigma}_\ell^2)}, and \eqn{\bar{\sigma}_\ell^2} {~} \eqn{Uniform(0, d)}, with \eqn{d} a hyper-parameter. \cr \eqn{\eta_i = \sum_{p = 0}^{P}\sum_{j=1}^{J_p}X_{i,j}^p\beta_{j,s_i}^p}, \eqn{s_i} is the subject id of response \eqn{i}. see \code{\link{BANOVA-package}}
  
}

\value{
  \code{BANOVA.ordMultinomial} returns an object of class \code{"BANOVA.ordMultinomial"}. The returned object is a list containing:
  
  \item{anova.table}{table of effect sizes \code{\link{BAnova}}}
  \item{coef.tables}{table of estimated coefficients}
  \item{pvalue.table}{table of p-values \code{\link{table.pvalues}}}
  \item{dMatrice}{design matrices at level 1 and level 2}
  \item{samples_l2_param}{posterior samples of level 2 parameters}
  \item{samples_cutp_param}{posterior samples of cutpoints}
  \item{data}{original data.frame}
  \item{mf1}{model.frame of level 1}
  \item{mf2}{model.frame of level 2}
  \item{JAGSmodel}{'JAGS' model}
}

\examples{
data(goalstudy)
\dontshow{
res <- BANOVA.ordMultinomial (perceivedsim~1, ~progress*prodvar,  goalstudy,
goalstudy$id, burnin = 5, sample = 5, thin = 1)
}
\donttest{
res <- BANOVA.ordMultinomial (perceivedsim~1, ~progress*prodvar,  goalstudy,
goalstudy$id, burnin = 5000, sample = 2000, thin = 10)
summary(res)
# or use BANOVA.run based on 'Stan'
res <- BANOVA.run(perceivedsim~1, ~progress*prodvar, data = goalstudy, 
model_name = 'ordMultinomial', id = goalstudy$id, iter = 100, thin = 1, chains = 2)
}

}
