\name{table.pvalues}
\alias{table.pvalues}
\title{Function to print the table of p-values}
\description{
Computes the Baysian p-values for the test concerning all coefficients/parameters: \cr

For \eqn{p = 1,...,P} \cr
\eqn{H_0:\theta_{j,k}^{p,q}=0} \cr
\eqn{H_1:\theta_{j,k}^{p,q} \neq 0} \cr

The two-sided P-value for the sample outcome is obtained by first finding the one sided P-value, \eqn{min(P(\theta_{j,k}^{p,q}<0),P(\theta_{j,k}^{p,q}>0 ))} which can be estimated from posterior samples. For example, \eqn{P(\theta_{j,k}^{p,q}>0) = \frac{n_+}{n}}, where \eqn{n_+} is the number of posterior samples that are greater than 0, \eqn{n} is the target sample size. The two sided P-value is \eqn{P_\theta(\theta_{j,k}^{p,q}) = 2*min(P(\theta_{j,k}^{p,q}<0),P(\theta_{j,k}^{p,q}>0 ))}.

If there are \eqn{\theta_{j,k_1}^{p,q},\theta_{j,k_2}^{p,q},...,\theta_{j,k_J}^{p,q}} representing J levels of a multi-level variable, we use a single P-value to represent the significance of all levels. The two alternatives are:

\eqn{H_0:\theta_{j,k_1}^{p,q} = \theta_{j,k_2}^{p,q} = \cdots = \theta_{j,k_J}^{p,q}=0} \cr
\eqn{H_1} : some \eqn{\theta_{j,k_j}^{p,q} \neq 0}

Let \eqn{\theta_{j,k_{min}}^{p,q}} and \eqn{\theta_{j,k_{max}}^{p,q}} denote the coefficients with the smallest and largest posterior mean. Then the overall P-value is defined as 

\eqn{min(P_\theta (\theta_{j,k_{min}}^{p,q}), P_\theta(\theta_{j,k_{max}}^{p,q}))}. 
}

\usage{
table.pvalues(x)
}
\arguments{
  \item{x}{the object from BANOVA.*}
}

\source{
It borrows the idea of Sheffe F-test for multiple testing: the F-stat for testing the contrast with maximal difference from zero. Thank Dr. P. Lenk of the University of Michigan for this suggestion.
}
\examples{

data(goalstudy)
res1 <- BANOVA.Normal(bid~1, ~progress*prodvar,  goalstudy, goalstudy$id, 
burnin = 1000, sample = 1000, thin = 2)
\donttest{
# or use BANOVA.run
res1 <- BANOVA.run(bid~1, ~progress*prodvar, model_name = "Normal", 
data = goalstudy, id = goalstudy$id, iter = 1000, thin = 1, chains = 2) 
}
table.pvalues(res1)

}

