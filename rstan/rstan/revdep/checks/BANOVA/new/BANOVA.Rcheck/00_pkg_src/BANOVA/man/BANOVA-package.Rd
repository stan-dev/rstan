\name{BANOVA-package}
\alias{BANOVA-package}
\alias{BANOVA}
\docType{package}
\title{
BANOVA: Hierarchical Bayesian ANOVA Models
}
\description{
This package includes several hierarchical Bayes Analysis of Variance models. These models are suited for the analysis of experimental designs in which both within- and between- subjects factors are manipulated, and account for a wide variety of distributions of the dependent variable. Floodlight analysis and mediation analysis basaed on these models are also provided. The package uses 'Stan' and 'JAGS' as the computational platform.
}
\details{
\tabular{ll}{
Package: \tab BANOVA\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 2019-01-21\cr
License: \tab GPL (>= 2)\cr
}
Model: \cr
  \eqn{E(y_i) = g^{-1}(\eta_i)} \cr
  where \eqn{\eta_i = \sum_{p = 0}^{P}\sum_{j=1}^{J_p}X_{i,j}^p\beta_{j,s_i}^p}, \eqn{s_i} is the subject id of data response \eqn{i}. Missing values (NAs) of \eqn{y_i} are allowed. The within-subjects factors and their interactions are indexed by \eqn{p (p = 1,2,.,P)}. Each index \eqn{p} represents a batch of \eqn{J_p} coefficients: \eqn{\beta_{j,s}^p, j = 1,.,J_p};\eqn{s = 1,.,S} indexes subjects. Note that if the subject-level covariate is continuous, \eqn{J_p=1}, so that ANCOVA models are also accommodated (relaxing their "constant slope" assumption).  \cr
  
  The population-level model allows for heterogeneity among subjects, because the subject-level coefficients \eqn{\beta_{j,s}^p} are assumed to follow a multivariate normal distribution.The between-subjects factors and their interactions are indexed by \eqn{q,(q = 1,2,.,Q)}, \eqn{q = 0} denotes the constant term. The population-level ANOVA can be written as: \cr
  \eqn{\beta_{j,s}^p = \sum_{q = 0}^Q \theta_{j,k_s^q}^{pq} + \delta_{j,s}^p} \cr
  
  The population-level ANCOVA model can be expressed as a linear model with a design matrix \eqn{Z} that contains all between-subjects factors and their interactions and a constant term: \cr
  \eqn{\beta_{j,s}^p = \sum_{k = 1}^Q Z_{s,k}\theta_{j,k}^{p} + \delta_{j,s}^p} \cr
  where \eqn{Z_{s,k} }is an element of \eqn{Z}, a \eqn{S \times Q} matrix of covariates. \eqn{\theta_{j,k}^p} is a hyperparameter which captures the effects of between-subjects factor \eqn{q} on the parameter \eqn{\beta_{j,s}^p} of within-subjects factor p. The error \eqn{\delta_{j,s}^p} is assumed to be normal: \eqn{\delta_{j,s}^p} {~} \eqn{N(0,\sigma_p^{-2} )}. Proper, but diffuse priors are assumed: \eqn{\theta_{j,k}^p} {~} \eqn{N(0,\gamma)}, and \eqn{\sigma_p^{-2}} {~} \eqn{Gamma(a,b)}, where \eqn{\gamma,a,b} are hyper-parameters. The default setting is \eqn{\gamma = 10^{-4}, a = 1, b = 1}. \cr
  
  Note that missing values of independent variables are currently not allowed in the package.
}
\author{
Chen Dong;
Michel Wedel

Maintainer: Chen Dong <cdong@math.umd.edu>

}
\references{
Dong, C. and Wedel, M. (2017) 
\emph{BANOVA: An R Package for Hierarchical Bayesian ANOVA}, Journal of Statistical Software, Vol. 81, No.9, pp. 1-46. 

McCullagh, P., Nelder, JA. (1989)
\emph{Generalized linear models}, New York, NY: Chapman and Hall.

Gelman, A. (2005)
\emph{Analysis of variance-why it is more important than ever}, Ann. Statist., Vol. 33, No. 1, pp. 1-53.

Rossi, P., Allenby,G., McCulloch, R. (2005)
\emph{Bayesian Statistics and Marketing}, John Wiley and Sons.

Gill, J. (2007) 
\emph{Bayesian Methods for the Social and Behavioral Sciences}, Chapman and Hall, Second Edition.

Gelman, A., Carlin, J., Stern, H. and Dunson, D. (2013)
\emph{Bayesian Data Analysis}, London: Chapman and Hall.

Wedel, M. and Dong, C. (2016) \emph{BANOVA: Bayesian Analysis of Variance for Consumer Research}. Submitted.

}
%Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
% directory ~~
\keyword{ Hierarchical Bayes, ANOVA }


