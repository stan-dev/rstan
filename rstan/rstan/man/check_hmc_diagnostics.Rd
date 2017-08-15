\name{check_hmc_diagnostics}
\alias{check_hmc_diagnostics}
\alias{check_divergences}
\alias{check_treedepth}
\alias{check_energy}

\title{Check HMC diagnostics after sampling} 

\description{
These functions print summaries of important HMC diagnostics.
The \code{check_hmc_diagnostics} function calls the other functions internally
and prints an overall summary, but the other functions can also be called 
directly: 
\code{check_divergences} prints the number (and percentage) of iterations that
ended with a divergence, \code{check_treedepth} prints the number (and percentage) of iterations that saturated the max treedepth, and \code{check_energy} 
prints E-BFMI values for each chain for which E-BFMI is less than 0.2.
} 
 
\usage{
check_hmc_diagnostics(object)
check_divergences(object)
check_treedepth(object)
check_energy(object)
} 
 
\arguments{
  \item{object}{A stanfit object.}
} 

\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{http://mc-stan.org/}. 
  
  Betancourt, M. (2017). A conceptual introduction to Hamiltonian Monte Carlo.
  \url{https://arxiv.org/abs/1701.02434}.
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
