\name{conv.diag}
\alias{conv.diag}
\title{Function to display the convergence diagnostics}
\description{
The Geweke diagnostic and the Heidelberg and Welch diagnostic are reported. These two convergence diagnostics are calculated based on only a single MCMC chain. Both diagnostics require a single chain and may be applied with any MCMC method. The functions \code{\link[coda]{geweke.diag}}, \code{\link[coda]{heidel.diag}} in \bold{coda} package is used to compute this diagnostic. 

Geweke's convergence diagnostic is calculated by taking the difference between the means from the first \eqn{n_A} iterations and the last \eqn{n_B} iterations. If the ratios \eqn{n_A/n}  and \eqn{n_B/n} are fixed and \eqn{nA +nB < n}, then by the central limit theorem, the distribution of this diagnostic approaches a standard normal as  \eqn{n} tends to infinity. In our package, \eqn{n_A= .2*n} and \eqn{n_B= .5*n}. 

The Heidelberg and Welch diagnostic is based on a test statistic to accept or reject the null hypothesis that the Markov chain is from a stationary distribution. The present package reports the stationary test.The convergence test uses the Cramer-von Mises statistic to test for stationary. The test is successively applied on the chain. If the null hypothesis is rejected, the first 10\% of the iterations are discarded and the stationarity test repeated. If the stationary test fails again, an additional 10\% of the iterations are discarded and the test repeated again. The process continues until 50\% of the iterations have been discarded and the test still rejects. In our package, \eqn{eps = 0.1, pvalue = 0.05} are used as parameters of the function \code{\link[coda]{heidel.diag}}.

}
\usage{
conv.diag(x)
}
\arguments{
  \item{x}{the object from BANOVA.*}

}
\value{
  \code{conv.diag} returns a list of two diagnostics:
  
  \item{sol_geweke}{The Geweke diagnostic}
  \item{sol_heidel}{The Heidelberg and Welch diagnostic}
}

\references{
Plummer, M., Best, N., Cowles, K. and Vines K. (2006) \emph{CODA: Convergence Diagnosis and Output Analysis for MCMC}, R News, Vol 6, pp. 7-11.

Geweke, \emph{J. Evaluating the accuracy of sampling-based approaches
to calculating posterior moments}, In \emph{Bayesian Statistics 4}
(ed JM Bernado, JO Berger, AP Dawid and AFM Smith). Clarendon Press,
Oxford, UK.

Heidelberger, P. and Welch, PD. (1981)
\emph{A spectral method for confidence interval generation and run length control in simulations}, Comm. ACM. Vol. 24, No.4, pp. 233-245. 

Heidelberger, P. and Welch, PD. (1983)
\emph{Simulation run length control in the 
presence of an initial transient}, Opns Res., Vol.31, No.6, pp. 1109-44.

Schruben, LW. (1982)
\emph{Detecting initialization bias in simulation experiments},
Opns. Res., Vol. 30, No.3, pp. 569-590.

}
\examples{
data(goalstudy)
\donttest{
res1 <- BANOVA.run(bid~1, ~progress*prodvar, model_name = "Normal", data = goalstudy, 
id = goalstudy$id, iter = 100, thin = 1) 
conv.diag(res1)
# might need pairs() to confirm the convergence
}
}