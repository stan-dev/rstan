
\name{stanfit-class} 
\Rdversion{1.1}

\docType{class}
\alias{stanfit-class}
\alias{stanfit}
\alias{summary,stanfit-method}
%% \alias{extract,stanfit-method}
\alias{show,stanfit-method}
% \alias{plot,stanfit-method}
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
\alias{get_posterior_mean}
\alias{get_posterior_mean,stanfit-method}
\alias{get_elapsed_time}
\alias{get_elapsed_time,stanfit-method}
\alias{get_inits,stanfit-method}
\alias{get_logposterior} 
\alias{get_logposterior,stanfit-method}
\alias{get_adaptation_info} 
\alias{get_adaptation_info,stanfit-method}
\alias{get_sampler_params}
\alias{get_sampler_params,stanfit,logical-method}
%% \alias{get_num_upars}
%% \alias{get_num_upars,stanfit-method}
%% \alias{unconstrain_pars}
%% \alias{unconstrain_pars,stanfit-method}
%%\alias{unconstrain_pars,stanfit,list-method}
%% \alias{constrain_pars}
%% \alias{constrain_pars,stanfit-method}
%%\alias{constrain_pars,stanfit,numeric-method}
%%\alias{log_prob}
%%\alias{log_prob,stanfit-method}
%% \alias{log_prob,stanfit,numeric-method}
%%\alias{grad_log_prob}
%%\alias{grad_log_prob,stanfit-method}
%% \alias{grad_log_prob,stanfit,numeric-method}
%% \alias{traceplot,stanfit-method}
\title{Class \code{stanfit}: fitted Stan model}

\description{
The output derived from fitting a Stan model, including the posterior draws,
as returned by the top-level function \code{\link{stan}} or the lower-level
methods \code{\link{sampling}} and \code{\link{vb}} defined on class
\code{stanmodel}. Methods \code{print}, \code{plot}, and many others, are
provided for summarizing results. Various access methods also allow the
underlying data (draws, diagnostics, etc.) comprising the fit to be retrieved.
Functions for computing the log probability density \code{lp__} and the
gradient are also exposed. The components (slots) of a \code{stanfit} object
and the various available methods are described below.
} 

\section{Objects from the Class}{
  Objects are created either by calling the \code{\link{stan}} function or
  one of the \code{\link{sampling}} and \code{\link{vb}} methods for S4 class 
  \code{\linkS4class{stanmodel}}.
}

\section{Slots}{
  \describe{
    \item{\code{model_name}:}{The model name as a string.} 
    \item{\code{model_pars}:}{A character vector of names of parameters
      (including transformed parameters and derived quantities).}
    \item{\code{par_dims}:}{A named list giving the dimensions for all 
      parameters. The dimension for a scalar parameter is given as 
      \code{numeric(0)}.}
    \item{\code{mode}:}{The mode of the fitted model, an \code{integer}. 
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
  \describe{
    \item{\code{show}}{\code{signature(object = "stanfit")}: print the default 
      summary for the model.}
    \item{\code{print}}{\code{signature(object = "stanfit")}: print a 
      customizable summary for the model. See the separate documentation 
      in \code{\link{print.stanfit}}.}
    \item{\code{plot}}{\code{signature(x = "stanfit", y = "missing")}: 
      create various plots summarizing the fitted model. See the separate 
      documentation on plotting at \code{\link{plot-methods}}.}
    \item{\code{summary}}{\code{signature(object = "stanfit")}: summarize the 
      distributions of all estimated quantities using the posterior draws.
      Included are quantiles 
      (e.g., 2.5\%, 50\%, 97.5\%, which can be specified by 
      using the optional \code{probs} argument in the same way as for 
      \code{\link[stats]{quantile}}), the mean, the standard deviation (sd), 
      the effective sample size (n_eff), and the split Rhat (the potential scale 
      reduction derived from all chains after splitting each chain in half and 
      treating the halves as chains).
      The \code{summary} method returns a named list with elements 
      \code{summary} and \code{c_summary}, which contain summaries for
      for all chains merged and individual chains, respectively.
      For the summary of all chains merged, also reported is \code{se_mean},
      the (Monte Carlo) standard error of the mean (not to be confused with the
      posterior standard deviation). 
      To specify parameters of interest, use the argument \code{pars}. The 
      default for \code{pars} is all parameters saved in the fitted
      results as well as the log-posterior (\code{lp__}).
      Another argument for \code{summary} is \code{use_cache}, which defaults
      to \code{TRUE}. When \code{use_cache=TRUE}, the summary quantities for
      all parameters are computed and cached for future use. 
      Setting \code{use_cache=FALSE} can be used to avoid performing the summary
      computations for all parameters if \code{pars} is given as some specific 
      parameters.}
    \item{\code{extract}}{\code{signature(object = "stanfit")}: 
      get the samples for all chains for all (or specified) parameters. 
      (see the separate documentation in \code{\link{extract}}). 
      Also see \code{\link{as.array.stanfit}} for coercing samples without 
      warmup to an array, matrix or data frame.}
    \item{\code{as.mcmc.list}}{\code{signature(object = "stanfit")}:
      return a list of all the chains that can be treated as an
      \code{\link[coda]{mcmc.list}} as in package \pkg{coda}. 
      Parameters of interest can be specified using \code{pars}. 
      The warmup samples are discarded. This method is not exported but
      \code{\link{As.mcmc.list}} is.}
    \item{\code{traceplot}}{\code{signature(object = "stanfit")}:
      plot the traces of the sample for all chains. See the documentation in 
      \code{\link{traceplot}}.} 
    \item{\code{get_posterior_mean}}{\code{signature(object = "stanfit")}:
      get the posterior mean for parameters of interest (using \code{pars}
      to specify) among \emph{all} parameters. Returned is a matrix with 
      one column per chain and an additional column for all chains combined.}
    \item{\code{get_elapsed_time}}{\code{signature(object = "stanfit")}:
      get the warmup time and sample time in seconds.
      A matrix of two columns is returned with each row containing the warmup
      and sample times for one chain.}
    \item{\code{get_stancode}}{\code{signature(object = "stanfit")}:
      get the Stan code for the fitted model as a string. The result can 
      be printed in a readable format using \code{\link{cat}}.}
    \item{\code{get_cppo_mode}}{\code{signature(object = "stanfit")}:
      get the optimization mode used for compilation. The returned string is 
      one of \code{"fast"}, \code{"presentation2"}, \code{"presentation1"},
      and \code{"debug"}.}
    \item{\code{get_stanmodel}}{\code{signature(object = "stanfit")}: 
      get the object of S4 class \code{\linkS4class{stanmodel}} of the fitted 
      model.} 
    \item{\code{get_inits}}{\code{signature(object = "stanfit")}:
      get the initial values for parameters used in sampling all chains. The 
      returned object is a list with the same structure as the \code{inits} 
      slot described above. If \code{object@mode=2} (error mode) an empty list 
      is returned.}
    \item{\code{get_seed}}{\code{signature(object = "stanfit")}: 
      get the (P)RNG seed used. When the fitted object 
      is empty (\code{mode=2}), \code{NULL} might be returned.
      In the case that the seeds for all chains are different, use 
      \code{get_seeds}.} 
    \item{\code{get_seeds}}{\code{signature(object = "stanfit")}: 
      get the seeds used for all chains. When the fitted object 
      is empty (\code{mode=2}), \code{NULL} might be returned.} 
    \item{\code{get_logposterior}}{\code{signature(object = "stanfit")}: 
      get the log-posterior (up to an additive constant, i.e. up to
      a multiplicative constant on the linear scale) for all chains. 
      Each element of the returned \code{list} is the log-posterior 
      for a single chain (at each iteration). 
      The optional argument \code{inc_warmup} (defaulting to \code{TRUE}) 
      indicates whether to include the warmup period.} 
    \item{\code{get_adaptation_info}}{\code{signature(object = "stanfit")}:
      obtain the adaptation information for the sampler if NUTS was used. 
      The results are returned as a list, each element of which is a character 
      string with the info for a single chain.} 
    \item{\code{get_sampler_params}}{\code{signature(object = "stanfit")}: 
      obtain the parameters used for the sampler such as 
      \code{stepsize} and \code{treedepth}. The results are returned
      as a list with one component (an array) per chain.
      The array has number of columns corresponding to the number
      of parameters used in the sampler and its column names 
      provide the parameter names. Optional argument \code{inc_warmup} 
      (defaulting to \code{TRUE}) indicates whether to include the warmup period.} 
    \item{\code{log_prob}}{\code{signature(object = "stanfit", "numeric")}: compute
      the log probability density (\code{lp__}) for a set of parameter 
      values (on the \emph{unconstrained} space) up to an additive constant. 
      The unconstrained parameters are specified using a numeric vector.
      The number of parameters on the unconstrained space can be obtained 
      using method \code{get_num_upars}. A numeric value is returned. See also
      the documentation in \code{\link{log_prob}}.} 
    \item{\code{grad_log_prob}}{\code{signature(object = "stanfit", "numeric")}: 
      compute the gradient of log probability density function for a set of parameter 
      values (on the \emph{unconstrained} space) up to an additive constant.
      The unconstrained parameters are specified using a numeric vector
      with the length being the number of unconstrained parameters.
      A numeric vector is returned with the length of the number of 
      unconstrained parameters and an attribute named \code{log_prob} being
      the \code{lp__}. See also the documentation in \code{\link{grad_log_prob}}.}
    \item{\code{get_num_upars}}{\code{signature(object = "stanfit")}: get the 
      number of unconstrained parameters of the model. The number of parameters 
      for a model is not necessarily equal to this number of unconstrained 
      parameters. For example, when a parameter is specified as a simplex of 
      length K, the number of unconstrained parameters is K-1.} 
    \item{\code{unconstrain_pars}}{\code{signature(object = "stanfit", "list")}: 
      transform the parameter to unconstrained space. The input is a named list 
      as for specifying initial values for each parameter. A numeric vector is
      returned. See also the documentation in \code{\link{unconstrain_pars}}.}
    \item{\code{constrain_pars}}{\code{signature(object = "stanfit", "numeric")}: 
      get the parameter values from their unconstrained space. The input is a 
      numeric vector. A list is returned.  This function is contrary to
      \code{unconstrain_pars}. See also the documentation in 
      \code{\link{constrain_pars}}.}
    \item{\code{dim}, \code{dimnames}, \code{names}}{See the documentation in 
      \code{\link{dim.stanfit}}.}
  }
}

\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{http://mc-stan.org}. 
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
plot(exfit)

adaptinfo <- get_adaptation_info(exfit)
seed <- get_seed(exfit)
sp <- get_sampler_params(exfit)
ml <- As.mcmc.list(exfit)
cat(get_stancode(exfit))
}

\keyword{classes}
