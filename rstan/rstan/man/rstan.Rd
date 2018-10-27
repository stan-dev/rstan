
\name{rstan-package}
\alias{rstan-package}
\alias{rstan} 
\docType{package}
\title{RStan --- the \R interface to Stan}
\description{
\if{html}{
   \figure{stanlogo.png}{options: width="50px" alt="mc-stan.org"}
   \emph{Stan Development Team}
}

  RStan is the \R interface to the \href{http://mc-stan.org/}{Stan} C++ package. 
  The RStan interface (\pkg{rstan} \R package) provides:
  \itemize{
  \item Full Bayesian inference using the No-U-Turn sampler (NUTS), 
    a variant of Hamiltonian Monte Carlo (HMC)
  \item Approximate Bayesian inference using automatic differentiation 
    variational inference (ADVI)
  \item Penalized maximum likelihood estimation using L-BFGS optimization
  }
  
  For documentation on Stan itself, including the manual and user guide
  for the modeling language, case studies and worked examples, and other
  tutorial information visit the Users section of the Stan website:
  
  \itemize{
  \item \href{http://mc-stan.org/users/documentation/}{mc-stan.org/users/documentation}
  }
  
}
\section{Other \R packages from the Stan Development Team}{
  Various related \R packages are also available from the Stan Development Team: 
  \tabular{llll}{
  \strong{Package} \tab \strong{Description} \tab \strong{Doc} \tab \strong{Website} \cr
  \pkg{bayesplot} 
    \tab ggplot-based plotting of parameter estimates, diagnostics, and posterior predictive checks. 
    \tab \link[bayesplot]{bayesplot-package}
    \tab \href{http://mc-stan.org/bayesplot/}{mc-stan.org/bayesplot}
    \cr
  \pkg{shinystan} 
    \tab Interactive GUI for exploring MCMC output. 
    \tab \link[shinystan]{shinystan-package}
    \tab \href{http://mc-stan.org/shinystan/}{mc-stan.org/shinystan}
    \cr
  \pkg{loo} 
    \tab Out-of-sample predictive performance estimates and model comparison.
    \tab \link[loo]{loo-package} 
    \tab \href{http://mc-stan.org/loo/}{mc-stan.org/loo}
    \cr
  \pkg{rstanarm} 
    \tab \R formula interface for applied regression modeling. 
    \tab \link[rstanarm]{rstanarm-package}
    \tab \href{http://mc-stan.org/rstanarm/}{mc-stan.org/rstanarm}
    \cr 
  \pkg{rstantools} 
    \tab Tools for developers of \R packages interfacing with Stan.
    \tab \link[rstantools]{rstantools-package}
    \tab \href{http://mc-stan.org/rstantools/}{mc-stan.org/rstantools}
    \cr
  }
}
\author{
  \tabular{ll}{
    Jonah Gabry (author)\tab <jonah.sol.gabry@columbia.edu> \cr
    Ben Goodrich (maintainer, author)\tab <benjamin.goodrich@columbia.edu> \cr
    Jiqiang Guo (author)\tab <guojq28@gmail.com> \cr 
  }

  There are also many other important contributors to RStan 
  (\href{https://github.com/stan-dev/rstan}{github.com/rstan}). 
  Please use 'Stan Development Team' whenever citing the R interface to Stan. 
  A BibTex entry is available from \url{http://mc-stan.org/rstan/authors} 
  or \code{citation("rstan")}.
}

\seealso{
\itemize{
  \item The RStan vignettes: \url{http://mc-stan.org/rstan/articles/}.
  \item \code{\link{stan}} for details on fitting models and 
  \code{\linkS4class{stanfit}} for information on the fitted model objects.
  \item The \code{\link{lookup}} for finding a function in the Stan language 
  that corresponds to a \R function or name.
  \item \url{https://github.com/stan-dev/rstan/issues/} to submit a bug
  report or feature request.
  \item \url{http://discourse.mc-stan.org} to ask a 
  question on the Stan Forums.
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

# extract samples 
e <- extract(fit, permuted = FALSE) # return a list of arrays 
mu <- e$mu 

arr <- as.array(fit) # return an array 
print(dimnames(arr))

mat <- as.matrix(fit) # return a matrix
print(dimnames(mat))
}}
\keyword{rstan} 
\keyword{package}
