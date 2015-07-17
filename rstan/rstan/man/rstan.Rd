
\name{rstan-package}
\alias{rstan-package}
\alias{rstan} 
\docType{package}
\title{RStan --- \R interface to Stan}
\description{
  \R interface to Stan, which is a C++ package for obtaining Bayesian
  inference using the No-U-turn sampler, a variant of Hamiltonian Monte
  Carlo (see \url{http://mc-stan.org/}).
}

\details{
  \tabular{ll}{
    Package: \tab rstan\cr
    Type: \tab Package\cr
    Version: \tab 2.7.0\cr
    Date: \tab Jul 17, 2015\cr
    License: \tab GPL-3 \cr
  }

  The RStan package provides an interface to Stan. 
  For more information on Stan and its modeling language,
  see the \emph{Stan Modeling Language User's Guide and Reference Manual},   
  which is available \url{http://mc-stan.org/}.

}
\author{
  \tabular{ll}{
  Author:\tab Jiqiang Guo <guojq28@gmail.com>   \cr 
         \tab Ben Goodrich <benjamin.goodrich@columbia.edu> \cr 
  } 

  \tabular{ll}{
  Maintainer:\tab Jiqiang Guo <benjamin.goodrich@columbia.edu>    \cr
  } 
}
\references{
  The Stan Development Team 
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
  mu ~ normal(0, 10);
  y ~ normal(mu, 1); 
} 
"

y <- rnorm(20) 
dat <- list(N = 20, y = y); 
fit <- stan(model_code = stanmodelcode, model_name = "example", 
            data = dat, iter = 2012, chains = 3, sample_file = 'norm.csv',
            verbose = TRUE) 
print(fit)

# extract samples 
e <- extract(fit, permuted = TRUE) # return a list of arrays 
mu <- e$mu 

m <- extract(fit, permuted = FALSE, inc_warmup = FALSE) # return an array 
print(dimnames(m))

# using as.array directly on stanfit objects 
m2 <- as.array(fit)

}}  

\seealso{
  \code{\link{stan}}, \code{\linkS4class{stanfit}} 
} 
\keyword{rstan} 
\keyword{package} 

