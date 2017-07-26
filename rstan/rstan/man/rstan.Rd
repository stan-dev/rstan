
\name{rstan-package}
\alias{rstan-package}
\alias{rstan} 
\docType{package}
\title{RStan --- \R interface to Stan}
\description{
\if{html}{
   \figure{stanlogo.png}{options: width="50px" alt="mc-stan.org"}
   \emph{Stan Development Team}
}

  RStan is the \R interface to the \href{http://mc-stan.org/}{Stan} C++ package. 
  RStan provides
  \itemize{
  \item full Bayesian inference using the No-U-Turn sampler (NUTS), 
    a variant of Hamiltonian Monte Carlo (HMC)
  \item approximate Bayesian inference using automatic differentiation 
    variational inference (ADVI)
  \item penalized maximum likelihood estimation using L-BFGS optimization
  }
  
  For documentation on the Stan modeling language see the
  \href{http://mc-stan.org/documentation/}{\emph{Stan Modeling Language User's Guide and Reference Manual}}.
}
\section{Other R packages from the Stan Development Team}{
  Various related \R packages are also available from the Stan Development Team: 
  \tabular{lll}{
  \strong{Package} \tab \strong{Description} \tab \strong{Link} \cr
  \pkg{bayesplot} 
    \tab ggplot-based plotting library for graphing parameter estimates, 
      MCMC diagnostics, and posterior predictive checks. 
    \tab \link[bayesplot]{bayesplot-package}
    \cr
  \pkg{shinystan} 
    \tab Interactive visual summaries and advanced posterior analysis of 
      MCMC output. 
    \tab \link[shinystan]{shinystan-package}
    \cr
  \pkg{loo} 
    \tab Out-of-sample predictive performance estimates and model comparison.
    \tab \link[loo]{loo-package} 
    \cr
  \pkg{rstanarm} 
    \tab \R formula interface for Bayesian applied regression modeling. 
    \tab \link[rstanarm]{rstanarm-package}
    \cr 
  \pkg{rstantools} 
    \tab Tools for developers of R packages interfacing with Stan.
    \tab \link[rstantools]{rstantools-package}
    \cr
  }
}
\author{
  \tabular{ll}{
  Authors:\tab Jiqiang Guo <guojq28@gmail.com>   \cr 
          \tab Ben Goodrich <benjamin.goodrich@columbia.edu> \cr
          \tab Jonah Gabry <jsg2201@columbia.edu> \cr
  }

  \tabular{ll}{
  Maintainer:\tab Ben Goodrich <benjamin.goodrich@columbia.edu>    \cr
  }
}

\seealso{
\itemize{
  \item \code{\link{stan}} for details on fitting models and 
  \code{\linkS4class{stanfit}} for information on the fitted model objects.
  \item \url{https://github.com/stan-dev/rstan/issues/} to submit a bug
  report or feature request.
  \item \url{http://discourse.mc-stan.org} to ask a 
  question on the Stan-users forum.
}
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
\keyword{rstan} 
\keyword{package}

