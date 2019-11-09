\name{BANOVA.Binomial}
\alias{BANOVA.Binomial}
\alias{predict.BANOVA.Binomial}
\alias{print.BANOVA.Binomial}
\alias{summary.BANOVA.Binomial}
\title{Estimation of BANOVA with a Binomial dependent variable}
\description{
\code{BANOVA.Binomial} implements a Hierarchical Bayesian ANOVA for a binomial response variable using a logit link and a normal heterogeneity distribution.

}
\usage{
BANOVA.Binomial(l1_formula = "NA", l2_formula = "NA", data,
  id, num_trials, l2_hyper = c(1, 1, 0.0001), burnin = 5000, sample = 2000, 
  thin = 10, adapt = 0, conv_speedup = F, jags = runjags.getOption('jagspath'))
\method{summary}{BANOVA.Binomial}(object, ...)
\method{predict}{BANOVA.Binomial}(object, newdata = NULL,...)
\method{print}{BANOVA.Binomial}(x, ...)
}
\arguments{
 \item{l1_formula}{formula for level 1 e.g. 'Y~X1+X2'}

  \item{l2_formula}{formula for level 2 e.g. '~Z1+Z2',
  response variable must not be included}

  \item{data}{a data.frame in long format including all features in level 1 and level 2(covariates and categorical factors) and responses}

  \item{id}{subject ID of each response unit}

  \item{num_trials}{the number of trials of each
  observation(=1, if it is bernoulli), the type is forced to be 'integer'}

  \item{l2_hyper}{level 2 hyperparameters, c(a, b, \eqn{\gamma}), default c(1,1,0.0001)}
  
  \item{burnin}{the number of burn in draws in the MCMC algorithm, default 5000}

  \item{sample}{target samples in the MCMC algorithm after thinning, default 2000}

  \item{thin}{the number of samples in the MCMC algorithm that needs to be thinned, default 10}
  
  \item{adapt}{the number of adaptive iterations, default 0 (see \link[runjags]{run.jags})}
  
  \item{conv_speedup}{whether to speedup convergence, default F}
  
  \item{jags}{the system call or path for activating 'JAGS'. Default calls findjags() to attempt to locate 'JAGS' on your system}
  
  \item{object}{object of class \code{BANOVA.Bin} (returned by \code{BANOVA.Bin})}

  \item{newdata}{test data, either a matrix, vector or a
  data frame. It must have the same format with the original data (the same column number)}
  
  \item{x}{object of class \code{BANOVA.Bin} (returned by \code{BANOVA.Bin})}
  
  \item{\dots}{additional arguments,currently ignored}
}

\details{
Level 1 model: \cr
  \eqn{y_i} {~} \eqn{Binomial(ntrials,p_i)}, \eqn{p_i = logit^{-1}(\eta_i)} \cr
  where ntrials is the binomial total for each record i, \eqn{\eta_i = \sum_{p = 0}^{P}\sum_{j=1}^{J_p}X_{i,j}^p\beta_{j,s_i}^p}, \eqn{s_i} is the subject id of response \eqn{i}. see \code{\link{BANOVA-package}}
  
}
\value{
  \code{BANOVA.Binomial} returns an object of class \code{"BANOVA.Bin"}. The returned object is a list containing:
  
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
data(colorad)
\donttest{
# mean center Blur for effect coding
colorad$blur <- colorad$blur - mean(colorad$blur)
res <- BANOVA.Binomial(y~typic, ~color*blur, colorad, colorad$id, as.integer(16), 
burnin = 5000, sample = 2000, thin = 10)
summary(res)
# or use BANOVA.run
res0 <- BANOVA.run(y~typic, ~color*blur, data = colorad, model_name = 'Binomial',
id = colorad$id, num_trials = as.integer(16), iter = 100, thin = 1, chains = 2)
summary(res0)
table.predictions(res)
# only in-model variables(except numeric variables) will be used
predict(res0, c(1, 0, 8, 2, 1, 0.03400759))
}

}
