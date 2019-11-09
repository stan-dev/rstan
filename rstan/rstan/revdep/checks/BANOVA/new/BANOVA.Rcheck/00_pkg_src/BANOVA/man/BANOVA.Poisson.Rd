\name{BANOVA.Poisson}
\alias{BANOVA.Poisson}
\alias{predict.BANOVA.Poisson}
\alias{print.BANOVA.Poisson}
\alias{summary.BANOVA.Poisson}
\title{Estimation of BANOVA with Poisson dependent variables}
\description{
\code{BANOVA.Poisson} implements a Hierarchical Bayesian ANOVA for models with a count-data response variable and normal heterogeneity distribution.
}
\usage{
BANOVA.Poisson(l1_formula = "NA", l2_formula = "NA",
  data, id, l2_hyper = c(1, 1, 0.0001), burnin = 5000, sample = 2000, thin = 10, 
  adapt = 0, conv_speedup = F, jags = runjags.getOption('jagspath'))
\method{summary}{BANOVA.Poisson}(object, ...)
\method{predict}{BANOVA.Poisson}(object, newdata = NULL,...)
\method{print}{BANOVA.Poisson}(x, ...)
}
\arguments{
  \item{l1_formula}{formula for level 1 e.g. 'Y~X1+X2'}

  \item{l2_formula}{formula for level 2 e.g. '~Z1+Z2',
  response variable must not be included, if missing, the single level model will be generated}

  \item{data}{a data.frame in long format including all features in level 1 and level 2(covariates and categorical factors) and responses}

  \item{id}{subject ID of each response unit}
  
  \item{l2_hyper}{level 2 hyperparameters, c(a, b, \eqn{\gamma}), default c(1,1,0.0001)}

  \item{burnin}{the number of burn in draws in the MCMC algorithm, default 5000}

  \item{sample}{target samples in the MCMC algorithm after thinning, default 2000}

  \item{thin}{the number of samples in the MCMC algorithm that needs to be thinned, default 10}
  
  \item{adapt}{the number of adaptive iterations, default 0 (see \link[runjags]{run.jags})}
  
  \item{conv_speedup}{whether to speedup convergence, default F}
  
  \item{jags}{the system call or path for activating 'JAGS'. Default calls findjags() to attempt to locate 'JAGS' on your system}
  
  \item{object}{object of class \code{BANOVA.Poisson} (returned by \code{BANOVA.Poisson})}

  \item{newdata}{test data, either a matrix, vector or a
  data frame. It must have the same format with the original data (the same column number)}
  
  \item{x}{object of class \code{BANOVA.Poisson} (returned by \code{BANOVA.Poisson})}
  
  \item{\dots}{additional arguments,currently ignored}
}
\details{
Level 1 model: \cr
  \eqn{y_i} {~} \eqn{Poisson(\lambda_i)}, \eqn{\lambda_i = exp(\eta_i + \epsilon_i)} \cr
  where \eqn{\eta_i = \sum_{p = 0}^{P}\sum_{j=1}^{J_p}X_{i,j}^p\beta_{j,s_i}^p}, \eqn{s_i} is the subject id of response \eqn{i}, see \code{\link{BANOVA-package}}. \eqn{\epsilon_i} is a dispersion term.
}

\value{
  \code{BANOVA.Poisson} returns an object of class \code{"BANOVA.Poisson"}. The returned object is a list containing:
  \item{anova.table}{table of effect sizes \code{\link{BAnova}}}
  \item{coef.tables}{table of estimated coefficients}
  \item{pvalue.table}{table of p-values \code{\link{table.pvalues}}}
  \item{dMatrice}{design matrices at level 1 and level 2}
  \item{samples_l2_param}{posterior samples of level 2 parameters}
  \item{samples_l2_sigma_param}{posterior samples of level 2 standard deviations}
  \item{data}{original data.frame}
  \item{mf1}{model.frame of level 1}
  \item{mf2}{model.frame of level 2}
  \item{JAGSmodel}{'JAGS' model}
}

\examples{
\donttest{
# use the bpndata dataset
data(bpndata)
# within-subjects model using the dependent variable : PIC_FIX 
res1 <- BANOVA.Poisson(PIC_FIX ~ AD_ID + PIC_SIZE+ PAGE_NUM 
+ PAGE_POS, ~1, bpndata, bpndata$RESPONDENT_ID, burnin = 500, 
sample = 200, thin = 5)
summary(res1)

# use the goalstudy dataset
data(goalstudy)
goalstudy$bid <- as.integer(goalstudy$bid + 0.5)
res2<-BANOVA.Poisson(bid~1, ~progress*prodvar,  goalstudy, goalstudy$id, 
burnin = 5000, sample = 2000, thin = 10)
summary(res2)

# or use the BANOVA.run based on 'Stan'
res3 <- BANOVA.run(bid~1, ~progress*prodvar, data = goalstudy, 
model_name = 'Poisson', id = goalstudy$id, iter = 100, thin = 1, chains = 2)

}
}


