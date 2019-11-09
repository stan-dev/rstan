\name{BANOVA.Bernoulli}
\alias{BANOVA.Bernoulli}
\alias{predict.BANOVA.Bernoulli}
\alias{print.BANOVA.Bernoulli}
\alias{summary.BANOVA.Bernoulli}
\title{Estimation of BANOVA with a Bernoulli dependent variable}
\description{
\code{BANOVA.Bernoulli} implements a Bayesian ANOVA for binary dependent variable, using a logit link and a normal heterogeneity distribution.
}
\usage{
BANOVA.Bernoulli(l1_formula = "NA", l2_formula = "NA", data,
  id, l2_hyper = c(1, 1, 0.0001), burnin = 5000, sample = 2000, thin = 10, 
  adapt = 0, conv_speedup = F, jags = runjags.getOption('jagspath'))

\method{summary}{BANOVA.Bernoulli}(object, ...)
\method{predict}{BANOVA.Bernoulli}(object, newdata = NULL,...)
\method{print}{BANOVA.Bernoulli}(x, ...)
}
\arguments{
  \item{l1_formula}{formula for level 1 e.g. 'Y~X1+X2'}

  \item{l2_formula}{formula for level 2 e.g. '~Z1+Z2',
  response variable must not be included}

  \item{data}{a data.frame in long format including all features in level 1 and level 2(covariates and categorical factors) and responses}

  \item{id}{subject ID of each response unit}
  
  \item{l2_hyper}{level 2 hyperparameters, c(a, b, \eqn{\gamma}), default c(1,1,0.0001)}

  \item{burnin}{the number of burn in draws in the MCMC algorithm, default 5000}

  \item{sample}{target samples in the MCMC algorithm after thinning, default 2000}

  \item{thin}{the number of samples in the MCMC algorithm that needs to be thinned, default 10}
  
  \item{adapt}{the number of adaptive iterations, default 0 (see \link[runjags]{run.jags})}
  
  \item{conv_speedup}{whether to speedup convergence, default F}
  
  \item{jags}{the system call or path for activating 'JAGS'. Default calls findjags() to attempt to locate 'JAGS' on your system}
  
  \item{object}{object of class \code{BANOVA.Bern} (returned by \code{BANOVA.Bern})}

  \item{newdata}{test data, either a matrix, vector or a
  data.frame. It must have the same format with the original data (the same number of features and the same data classes)}
  
  \item{x}{object of class \code{BANOVA.Bern} (returned by \code{BANOVA.Bern})}
  
  \item{\dots}{additional arguments,currently ignored}
}

\details{
Level 1 model: \cr
  \eqn{y_i} {~} \eqn{Binomial(1,p_i)}, \eqn{p_i = logit^{-1}(\eta_i)} \cr
  where \eqn{\eta_i = \sum_{p = 0}^{P}\sum_{j=1}^{J_p}X_{i,j}^p\beta_{j,s_i}^p}, \eqn{s_i} is the subject id of data record \eqn{i}. see \code{\link{BANOVA-package}}
  
}
\value{
  \code{BANOVA.Bernoulli} returns an object of class \code{"BANOVA.Bernoulli"}. The returned object is a list containing:
  
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
data(bernlogtime)
\donttest{
# model with the dependent variable : response
res <- BANOVA.Bernoulli(response~typical, ~blur + color, bernlogtime,
bernlogtime$subject, burnin = 5000, sample = 2000, thin = 10)
summary(res)
# or use BANOVA.run
res0 <- BANOVA.run(response~typical, ~blur + color, data = bernlogtime, 
model_name = 'Bernoulli',id = bernlogtime$subject, iter = 100, thin = 1, chains = 2)
summary(res0)
}

}

