\name{BANOVA.Multinomial}
\alias{BANOVA.Multinomial}
\alias{predict.BANOVA.Multinomial}
\alias{print.BANOVA.Multinomial}
\alias{summary.BANOVA.Multinomial}
\title{Estimation of BANOVA with a Multinomial dependent variable}
\description{
\code{BANOVA.Multinomial} implements a Hierarchical Bayesian ANOVA for multinomial response variable using a logit link and a normal heterogeneity distribution.
}
\usage{
BANOVA.Multinomial(l1_formula = "NA", l2_formula = "NA",
  dataX, dataZ, y, id, l2_hyper = c(1, 1, 0.0001), burnin = 5000, sample = 2000, 
  thin = 10, adapt = 0, conv_speedup = F, jags = runjags.getOption('jagspath'))
\method{summary}{BANOVA.Multinomial}(object, ...)
\method{predict}{BANOVA.Multinomial}(object, Xsamples = NULL, Zsamples = NULL,...)
\method{print}{BANOVA.Multinomial}(x, ...)
}
\arguments{
  \item{l1_formula}{formula for level 1 e.g. '~X1+X2',
  response variable must not be included}

  \item{l2_formula}{formula for level 2 e.g. '~Z1+Z2',
  response variable must not be included}

  \item{dataX}{a list of data frames(each corresponds to the choice set of each observation) that includes all covariates and factors}
  
  \item{dataZ}{a data frame(long format) that includes all level 2 covariates and factors}
  
  \item{y}{choice responses, 1,2,3...}

  \item{id}{subject id}

  \item{l2_hyper}{level 2 hyperparameters, c(a, b, \eqn{\gamma}), default c(1,1,0.0001)}
  
  \item{burnin}{the number of burn in draws in the MCMC algorithm, default 5000}

  \item{sample}{target samples in the MCMC algorithm after thinning, default 2000}

  \item{thin}{the number of samples in the MCMC algorithm that needs to be thinned, default 10}
  
  \item{adapt}{the number of adaptive iterations, default 0 (see \link[runjags]{run.jags})}
  
  \item{conv_speedup}{whether to speedup convergence, default F}
  
  \item{jags}{the system call or path for activating 'JAGS'. Default calls findjags() to attempt to locate 'JAGS' on your system}
  
  \item{object}{object of class \code{BANOVA.Multinomial}(returned by \code{BANOVA.Multinomial}) }

  \item{Xsamples}{new data samples in level one, must be a
  list( the same format with the traning data), numeric
  variables must be mean centered.}

  \item{Zsamples}{new data samples in level two( the same
  format with the traning data), numeric variables must be
  mean centered.}
  
  \item{x}{object of class \code{BANOVA.Multinomial} (returned by \code{BANOVA.Multinomial})}
  
  \item{\dots}{additional arguments,currently ignored}
}
\details{
Level 1 model: \cr
\eqn{P(y_i = \ell) = \frac{exp(\eta_{i\ell})}{\sum_{\ell=1}^L exp(\eta_{i\ell})}} \cr
where \eqn{\eta_{i\ell} = \sum_{p = 0}^{P}\sum_{j=1}^{J_p}X_{i,j}^{k,p} \beta_{j,s_i}^p}, \eqn{s_i} is the subject id of response \eqn{i}, see \code{\link{BANOVA-package}}. \eqn{X_{i,j}^{k,p}} is the design matrix corresponding to each class \eqn{\ell(\ell=1,.,L)} of \eqn{y_i}. 
The first level of the response is the base level, thus the intercept corresponding to this level will not be included.
}

\value{
  \code{BANOVA.Multinomial} returns an object of class \code{"BANOVA.Multinomial"}. The returned object is a list containing:
  
  \item{anova.table}{table of effect sizes \code{\link{BAnova}}}
  \item{coef.tables}{table of estimated coefficients}
  \item{pvalue.table}{table of p-values \code{\link{table.pvalues}}}
  \item{dMatrice}{design matrices at level 1 and level 2}
  \item{samples_l2_param}{posterior samples of level 2 parameters}
  \item{dataX}{original dataX}
  \item{dataZ}{original dataZ}
  \item{mf1}{model.frame of level 1}
  \item{mf2}{model.frame of level 2}
  \item{n_categories}{the number of categories of the response}
  \item{JAGSmodel}{'JAGS' model}
}

\examples{

# see 'choicedata'
data(choicedata)
# generate dataX(convert the within-subject variables to a list)
dataX <- list()
for (i in 1:nrow(choicedata)){
  logP <- as.numeric(log(choicedata[i,3:8]))
  # all numeric variables must be mean centered
  dataX[[i]] <- as.data.frame(logP) - mean(logP)
}
dataZ <- choicedata[,9:13]
\donttest{
res <- BANOVA.Multinomial(~ logP, ~ college, dataX, dataZ, 
  choicedata$choice, choicedata$hhid, burnin = 5000, sample = 2000, thin = 10)
# or use BANOVA.run based on 'Stan'
res <- BANOVA.run(~ logP, ~ college, dataX = dataX, dataZ = dataZ, 
                    model_name = 'Multinomial', y_value = choicedata$choice, 
                    id = choicedata$hhid, iter = 100, thin = 1, chains = 2)
}
}