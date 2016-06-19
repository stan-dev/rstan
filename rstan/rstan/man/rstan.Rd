
\name{rstan-package}
\alias{rstan-package}
\alias{rstan} 
\docType{package}
\title{RStan --- \R interface to Stan}
\description{
  RStan is the \R interface to the Stan C++ package. RStan provides
  \itemize{
  \item full Bayesian inference using the No-U-Turn sampler (NUTS), 
    a variant of Hamiltonian Monte Carlo (HMC)
  \item approximate Bayesian inference using automatic differentiation 
    variational inference (ADVI)
  \item penalized maximum likelihood estimation using L-BFGS optimization
  }
  
  For more information about Stan visit \url{http://mc-stan.org/}.
}

\details{
  \tabular{ll}{
    Package: \tab rstan\cr
    Version: \tab 2.10.0\cr
    Date: \tab June 20, 2016\cr
    License: \tab GPL-3 \cr
  }
  
  For more information on Stan and its modeling language,
  see the \emph{Stan Modeling Language User's Guide and Reference Manual}   
  available at \url{http://mc-stan.org/}.

}
\author{
  \tabular{ll}{
  Authors:\tab Jiqiang Guo <guojq28@gmail.com>   \cr 
          \tab Ben Goodrich <benjamin.goodrich@columbia.edu> \cr
          \tab Jonah Gabry >jsg2201@columbia.edu> \cr
  }

  \tabular{ll}{
  Maintainer:\tab Ben Goodrich <benjamin.goodrich@columbia.edu>    \cr
  } 
}
\references{
  Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{http://mc-stan.org/}. 
}

\examples{\dontrun{ 

stanmodelcode <- "
data {
  int<lower=0> N;
  real y[N];
} 

parameters {
  real mu;
} 

model {
  target += normal_lpdf(mu | 0, 10);
  target += normal_lpdf(y  | mu, 1);
} 
"

y <- rnorm(20) 
dat <- list(N = 20, y = y); 
fit <- stan(model_code = stanmodelcode, model_name = "example", 
            data = dat, iter = 2012, chains = 3, sample_file = 'norm.csv',
            verbose = TRUE) 
print(fit)
traceplot(fit)

# extract samples 
e <- extract(fit, permuted = TRUE) # return a list of arrays 
mu <- e$mu 

m <- extract(fit, permuted = FALSE, inc_warmup = FALSE) # return an array 
print(dimnames(m))

# using as.array directly on stanfit objects 
m2 <- as.array(fit)

}}  

\seealso{
  The \code{\link{stan}} function for details on fitting models and 
  \code{\linkS4class{stanfit}} for information on the fitted model objects.
  
  Several related \R packages are also available from the Stan Development Team: 
  \pkg{loo} (\link[loo]{loo-package}) offers model comparison on estimated
  out-of-sample predictive performance, 
  \pkg{shinystan} (\link[shinystan]{shinystan-package}) provides the ShinyStan 
  GUI for exploring fitted Bayesian models, 
  and \pkg{rstanarm} is an appendage to \pkg{rstan} providing an \R formula 
  interface for Bayesian regression modeling.
}
\keyword{rstan} 
\keyword{package} 

