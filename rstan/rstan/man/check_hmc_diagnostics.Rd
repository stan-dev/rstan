\name{check_hmc_diagnostics}
\alias{check_hmc_diagnostics}
\alias{check_divergences}
\alias{check_treedepth}
\alias{check_energy}
\alias{get_divergent_iterations}
\alias{get_max_treedepth_iterations}
\alias{get_num_leapfrog_per_iteration}
\alias{get_num_divergent}
\alias{get_num_max_treedepth}
\alias{get_bfmi}
\alias{get_low_bfmi_chains}


\title{Check HMC diagnostics after sampling} 

\description{
These functions print summaries of important HMC diagnostics or extract
those diagnostics from a \code{stanfit} object. See the \strong{Details}
section, below.
}

\usage{
check_hmc_diagnostics(object)
check_divergences(object)
check_treedepth(object)
check_energy(object)

get_divergent_iterations(object)
get_max_treedepth_iterations(object)
get_num_leapfrog_per_iteration(object)

get_num_divergent(object)
get_num_max_treedepth(object)

get_bfmi(object)
get_low_bfmi_chains(object)
} 
 
\arguments{
  \item{object}{A stanfit object.}
} 

\details{
The \code{check_hmc_diagnostics} function calls the other \code{check_*} 
functions internally and prints an overall summary, but the other 
functions can also be called directly: 

\itemize{
\item \code{check_divergences} prints the number (and percentage) of 
iterations that ended with a divergence,
\item \code{check_treedepth} prints the number (and percentage) of iterations 
that saturated the max treedepth, 
\item \code{check_energy} prints E-BFMI values for each chain for which E-BFMI 
is less than 0.2. 
}

The \code{get_*} functions are for programmatic access to the diagnostics. 

\itemize{
\item \code{get_divergent_iterations} and \code{get_max_treedepth_iterations} 
return a logical vector indicating problems for individual iterations,
\item \code{get_num_divergent} and \code{get_num_max_treedepth} return the 
number of offending interations, 
\item \code{get_num_leapfrog_per_iteration} returns an integer vector with the 
number of leapfrog evalutions for each iteration,
\item \code{get_bfmi} returns per-chain E-BFMI values and \code{get_low_bfmi_chains}
returns the indices of chains with low E-BFMI.
}

The following are several of many resources that provide more information on 
these diagnostics:
\itemize{
\item Brief explanations of some of the problems detected by these diagnostics 
can be found in the \href{https://mc-stan.org/misc/warnings.html}{\emph{Brief Guide to Stan's Warnings}}.

\item Betancourt (2017) provides much more depth on these diagnostics
as well as a conceptual introduction to Hamiltonian Monte Carlo in general.

\item Gabry et al. (2018) and the \pkg{bayesplot}
 package \href{https://mc-stan.org/bayesplot/articles/}{vignettes} demonstrate
 various visualizations of these diagnostics that can be made in \R.
}
}

\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{https://mc-stan.org/}. 
  
  Betancourt, M. (2017). A conceptual introduction to Hamiltonian Monte Carlo.
  \url{https://arxiv.org/abs/1701.02434}.
  
  Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A. (2018). 
  Visualization in Bayesian workflow. 
  \emph{Journal of the Royal Statistical Society Series A}, accepted for publication.
  arXiv preprint: https://arxiv.org/abs/1709.01449.
} 

\examples{
\dontrun{
schools <- stan_demo("eight_schools")
check_hmc_diagnostics(schools)
check_divergences(schools)
check_treedepth(schools)
check_energy(schools)
}
} 
