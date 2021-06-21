
\name{stanfit-class} 
\Rdversion{1.1}

\docType{class}
\alias{stanfit-class}
\alias{stanfit}
\alias{show,stanfit-method}
\alias{get_cppo_mode} 
\alias{get_cppo_mode,stanfit-method} 
\alias{get_stancode}
\alias{get_stancode,stanfit-method}
\alias{get_stanmodel}
\alias{get_stanmodel,stanfit-method}
\alias{get_seed}
\alias{get_seed,stanfit-method}
\alias{get_seeds}
\alias{get_seeds,stanfit-method}
\alias{get_inits}
\alias{get_inits,stanfit-method}
\alias{get_posterior_mean}
\alias{get_posterior_mean,stanfit-method}
\alias{get_elapsed_time}
\alias{get_elapsed_time,stanfit-method}
\alias{get_logposterior} 
\alias{get_logposterior,stanfit-method}
\alias{get_adaptation_info} 
\alias{get_adaptation_info,stanfit-method}
\alias{get_sampler_params}
\alias{get_sampler_params,stanfit,logical-method}

\title{Class \code{stanfit}: fitted Stan model}

\description{
The components (slots) of a \code{stanfit} object and the various available
methods are described below. When methods have their own more detailed 
documentation pages links are provided. 
} 

\section{Objects from the Class}{
An object of class \code{stanfit} contains the
output derived from fitting a Stan model as returned by the top-level function
\code{\link{stan}} or the lower-level methods \code{\link{sampling}} and
\code{\link{vb}} (which are defined on class \code{\linkS4class{stanmodel}}).
Many methods (e.g., \code{print}, \code{plot}, \code{summary}) are provided for 
summarizing results and various access methods also allow the underlying data 
(e.g., simulations, diagnostics) contained in the object to be retrieved.
}

\section{Slots}{
  \describe{
    \item{\code{model_name}:}{The model name as a string.} 
    \item{\code{model_pars}:}{A character vector of names of parameters
      (including transformed parameters and derived quantities).}
    \item{\code{par_dims}:}{A named list giving the dimensions for all 
      parameters. The dimension for a scalar parameter is given as 
      \code{numeric(0)}.}
    \item{\code{mode}:}{An integer indicating the mode of the fitted model.
       \code{0} indicates sampling mode, \code{1} indicates test gradient mode 
       (no sampling is done), and \code{2} indicates error mode (an error occurred
       before sampling). Most methods for \code{stanfit} objects are useful only 
       if \code{mode=0}.}
    \item{\code{sim}:}{A list containing simulation results including the 
      posterior draws as well as various pieces of metadata used by many of the 
      methods for \code{stanfit} objects.} 
    \item{\code{inits}:}{The initial values (either user-specified or generated
      randomly) for all chains. This is a list with one component per chain. Each
      component is a named list containing the initial values for each parameter
      for the corresponding chain.}
    \item{\code{stan_args}:}{A list with one component per chain containing the 
      arguments used for sampling (e.g. \code{iter}, \code{seed}, etc.).}
    \item{\code{stanmodel}:}{The instance of S4 class \code{\linkS4class{stanmodel}}.} 
    \item{\code{date}:}{A string containing the date and time the object was created.}
    \item{\code{.MISC}:}{Miscellaneous helper information used for the fitted model. 
      This is an object of type \code{environment}. Users rarely (if ever)
      need to access the contents of \code{.MISC}.}
  }
}

\section{Methods}{
  \strong{Printing, plotting, and summarizing:}
  \describe{
    \item{\code{show}}{Print the default summary for the model.}
    \item{\code{print}}{Print a customizable summary for the model. 
      See \code{\link{print.stanfit}}.}
    \item{\code{plot}}{Create various plots summarizing the fitted model. 
      See \code{\link{plot,stanfit-method}}.}
    \item{\code{summary}}{Summarize the distributions of estimated 
      parameters and derived quantities using the posterior draws. 
      See \code{\link{summary,stanfit-method}}.}
    \item{\code{get_posterior_mean}}{
      Get the posterior mean for parameters of interest (using \code{pars}
      to specify a subset of parameters). Returned is a matrix with 
      one column per chain and an additional column for all chains combined.}
  }
  \strong{Extracting posterior draws:}
  \describe{
    \item{\code{extract}}{Extract the draws for all chains for all 
      (or specified) parameters. See \code{\link{extract}}.}
    \item{\code{as.array}, \code{as.matrix}, \code{as.data.frame}}{
      Coerce the draws (without warmup) to an array, 
      matrix or data frame. See \code{\link{as.array.stanfit}}.}
    \item{\code{As.mcmc.list}}{Convert a \code{stanfit} object to an
      \code{\link[coda]{mcmc.list}} as in package \pkg{coda}. 
      See \code{\link{As.mcmc.list}}.}
    \item{\code{get_logposterior}}{
      Get the log-posterior at each iteration.
      Each element of the returned \code{list} is the vector of log-posterior
      values (up to an additive constant, i.e. up to a multiplicative constant 
      on the linear scale) for a single chain.
      The optional argument \code{inc_warmup} (defaulting to \code{TRUE}) 
      indicates whether to include the warmup period.}
  }
  \strong{Diagnostics, log probability, and gradients:}
  \describe{
    \item{\code{get_sampler_params}}{
      Obtain the parameters used for the sampler such as 
      \code{stepsize} and \code{treedepth}. The results are returned
      as a list with one component (an array) per chain.
      The array has number of columns corresponding to the number
      of parameters used in the sampler and its column names 
      provide the parameter names. Optional argument \code{inc_warmup} 
      (defaulting to \code{TRUE}) indicates whether to include the warmup period.}
    \item{\code{get_adaptation_info}}{
      Obtain the adaptation information for the sampler if NUTS was used. 
      The results are returned as a list, each element of which is a character 
      string with the info for a single chain.} 
    \item{\code{log_prob}}{
      Compute the log probability density (\code{lp__}) for a set of parameter 
      values (on the \emph{unconstrained} space) up to an additive constant. 
      The unconstrained parameters are specified using a numeric vector.
      The number of parameters on the unconstrained space can be obtained 
      using method \code{get_num_upars}. A numeric value is returned. See also
      the documentation in \code{\link{log_prob}}.} 
    \item{\code{grad_log_prob}}{
      Compute the gradient of log probability density function for a set of parameter 
      values (on the \emph{unconstrained} space) up to an additive constant.
      The unconstrained parameters are specified using a numeric vector
      with the length being the number of unconstrained parameters.
      A numeric vector is returned with the length of the number of 
      unconstrained parameters and an attribute named \code{log_prob} being
      the \code{lp__}. See also the documentation in \code{\link{grad_log_prob}}.}
    \item{\code{get_num_upars}}{
      Get the number of unconstrained parameters of the model. The number of parameters 
      for a model is not necessarily equal to this number of unconstrained 
      parameters. For example, when a parameter is specified as a simplex of 
      length K, the number of unconstrained parameters is K-1.} 
    \item{\code{unconstrain_pars}}{
      Transform the parameters to unconstrained space. The input is a named list 
      as for specifying initial values for each parameter. A numeric vector is
      returned. See also the documentation in \code{\link{unconstrain_pars}}.}
    \item{\code{constrain_pars}}{
      Get the parameter values from their unconstrained space. The input is a 
      numeric vector. A list is returned.  This function is contrary to
      \code{unconstrain_pars}. See also the documentation in 
      \code{\link{constrain_pars}}.}
  }
  \strong{Metadata and miscellaneous:}
  \describe{
    \item{\code{get_stancode}}{
      Get the Stan code for the fitted model as a string. The result can 
      be printed in a readable format using \code{\link{cat}}.}
    \item{\code{get_stanmodel}}{
      Get the object of S4 class \code{\linkS4class{stanmodel}} of the fitted 
      model.}
    \item{\code{get_elapsed_time}}{
      Get the warmup time and sample time in seconds.
      A matrix of two columns is returned with each row containing the warmup
      and sample times for one chain.}
    \item{\code{get_inits, iter = NULL}}{
      Get the initial values for parameters used in sampling all chains. The 
      returned object is a list with the same structure as the \code{inits} 
      slot described above. If \code{object@mode=2} (error mode) an empty list 
      is returned. If \code{iter} is not \code{NULL}, then the draw from that
      iteration is returned for each chain rather than the initial state.}
    \item{\code{get_cppo_mode}}{
      Get the optimization mode used for compilation. The returned string is 
      one of \code{"fast"}, \code{"presentation2"}, \code{"presentation1"},
      and \code{"debug"}.}
    \item{\code{get_seed}}{
      Get the (P)RNG seed used. When the fitted object 
      is empty (\code{mode=2}), \code{NULL} might be returned.
      In the case that the seeds for all chains are different, use 
      \code{get_seeds}.} 
    \item{\code{get_seeds}}{
      Get the seeds used for all chains. When the fitted object 
      is empty (\code{mode=2}), \code{NULL} might be returned.} 
  }
}

\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{https://mc-stan.org}.
}
\seealso{
  \code{\link{stan}} and \code{\linkS4class{stanmodel}} 
}
\examples{\dontrun{
showClass("stanfit")
ecode <- '
  parameters {
    real<lower=0> y[2];
  } 
  model {
    y ~ exponential(1);
  }
'
fit <- stan(model_code = ecode, iter = 10, chains = 1)
fit2 <- stan(fit = fit)
print(fit2)
plot(fit2)
traceplot(fit2)
ainfo <- get_adaptation_info(fit2)
cat(ainfo[[1]])
seed <- get_seed(fit2)
sp <- get_sampler_params(fit2)
sp2 <- get_sampler_params(fit2, inc_warmup = FALSE)
head(sp[[1]])

lp <- log_prob(fit, c(1, 2))
grad <- grad_log_prob(fit, c(1, 2))
lp2 <- attr(grad, "log_prob") # should be the same as "lp"

# get the number of parameters on the unconstrained space
n <- get_num_upars(fit)

# parameters on the positive real line (constrained space) 
y1 <- list(y = rep(1, 2)) 

uy <- unconstrain_pars(fit, y1) 
## uy should be c(0, 0) since here the log transformation is used
y1star <- constrain_pars(fit, uy)

print(y1)
print(y1star) # y1start should equal to y1 
}

# Create a stanfit object from reading CSV files of samples (saved in rstan
# package) generated by funtion stan for demonstration purpose from model as follows. 
# 
excode <- '
  transformed data {
    real y[20];
    y[1] <- 0.5796;  y[2]  <- 0.2276;   y[3] <- -0.2959; 
    y[4] <- -0.3742; y[5]  <- 0.3885;   y[6] <- -2.1585;
    y[7] <- 0.7111;  y[8]  <- 1.4424;   y[9] <- 2.5430; 
    y[10] <- 0.3746; y[11] <- 0.4773;   y[12] <- 0.1803; 
    y[13] <- 0.5215; y[14] <- -1.6044;  y[15] <- -0.6703; 
    y[16] <- 0.9459; y[17] <- -0.382;   y[18] <- 0.7619;
    y[19] <- 0.1006; y[20] <- -1.7461;
  }
  parameters {
    real mu;
    real<lower=0, upper=10> sigma;
    vector[2] z[3];
    real<lower=0> alpha;
  } 
  model {
    y ~ normal(mu, sigma);
    for (i in 1:3) 
      z[i] ~ normal(0, 1);
    alpha ~ exponential(2);
  } 
'

# exfit <- stan(model_code = excode, save_dso = FALSE, iter = 200, 
#               sample_file = "rstan_doc_ex.csv")
# 

exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'),
                       pattern='rstan_doc_ex_[[:digit:]].csv',
                       full.names = TRUE))

print(exfit)
\dontrun{
plot(exfit)
}

adaptinfo <- get_adaptation_info(exfit)
inits <- get_inits(exfit) # empty
inits <- get_inits(exfit, iter = 101)
seed <- get_seed(exfit)
sp <- get_sampler_params(exfit)
ml <- As.mcmc.list(exfit)
cat(get_stancode(exfit))
}

\keyword{classes}
