\name{stan}
\alias{stan}
\title{Fit a model with Stan}
\description{
\if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
Fit a model defined in the Stan modeling language and  
return the fitted result as an instance of \code{stanfit}. 
}
\usage{
stan(file, model_name = "anon_model", model_code = "", fit = NA, 
  data = list(), pars = NA,
  chains = 4, iter = 2000, warmup = floor(iter/2), thin = 1, 
  init = "random", seed = sample.int(.Machine$integer.max, 1), 
  algorithm = c("NUTS", "HMC", "Fixed_param"), %, "Metropolis"),
  control = NULL, sample_file = NULL, diagnostic_file = NULL, 
  save_dso = TRUE, verbose = FALSE, include = TRUE,
  cores = getOption("mc.cores", 1L),
  open_progress = interactive() && !isatty(stdout()) &&
                  !identical(Sys.getenv("RSTUDIO"), "1"),
  ...,
  boost_lib = NULL, eigen_lib = NULL
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
    The path to the Stan program to use. 
    \code{file} should be a character string file name or a connection that 
    \R supports containing the text of a model specification in the Stan 
    modeling language. 
    
    A model may also be specified directly as a character string using the 
    \code{model_code} argument, but we recommend always putting Stan programs 
    in separate files with a \code{.stan} extension.
    
    The \code{stan} function can also use the Stan program from an existing
    \code{stanfit} object via the \code{fit} argument. When \code{fit} is 
    specified, the \code{file} argument is ignored.
    }
    
  \item{model_code}{
    A character string either containing the model definition or the name of 
    a character string object in the workspace. This argument is used only 
    if arguments \code{file} and \code{fit} are not specified.
  }
    
  \item{fit}{
    An instance of S4 class \code{stanfit} derived from a previous fit;
    defaults to \code{NA}. If \code{fit} is not \code{NA}, the compiled model 
    associated with the fitted result is re-used; thus the time that would 
    otherwise be spent recompiling the C++ code for the model can be saved.
  }
    
  \item{model_name}{
    A string to use as the name of the model; defaults 
    to \code{"anon_model"}. However, the model name will be derived from 
    \code{file} or \code{model_code} (if \code{model_code} is the name
    of a character string object) if \code{model_name} is not specified.
    This is not a particularly important argument, although since it 
    affects the name used in printed messages, developers of other packages 
    that use \pkg{rstan} to fit models may want to use informative names.
  }
    
  \item{data}{A named \code{list} or \code{environment} providing the data for 
    the model, or a character vector for all the names of objects to use as data. 
    See the \strong{Passing data to Stan} section below.}
    
  \item{pars}{A character vector specifying parameters of interest to be saved. 
    The default is to save all parameters from the model. 
    If \code{include = TRUE}, only samples for parameters named in \code{pars} 
    are stored in the fitted results. Conversely, if \code{include = FALSE}, 
    samples for all parameters \emph{except} those named in \code{pars} are 
    stored in the fitted results.}
    
  \item{include}{Logical scalar defaulting to \code{TRUE} indicating
    whether to include or exclude the parameters given by the 
    \code{pars} argument. If \code{FALSE}, only entire multidimensional
    parameters can be excluded, rather than particular elements of them.}
    
  \item{iter}{A positive integer specifying the number of iterations for each 
    chain (including warmup). The default is 2000.}

  \item{warmup}{A positive integer specifying the number of warmup (aka burnin)
    iterations per chain. If step-size adaptation is on (which it is by default), 
    this also controls the number of iterations for which adaptation is run (and
    hence these warmup samples should not be used for inference). The number of 
    warmup iterations should be smaller than \code{iter} and the default is
    \code{iter/2}.}
    
  \item{chains}{A positive integer specifying the number of Markov chains. 
    The default is 4.} 
    
  \item{cores}{
    The number of cores to use when executing the Markov chains in parallel.
    The default is to use the value of the \code{"mc.cores"} option if it 
    has been set and otherwise to default to 1 core. However, we recommend 
    setting it to be as many processors as the hardware and RAM allow 
    (up to the number of chains). See \code{\link[parallel]{detectCores}} 
    if you don't know this number for your system.
  }

  \item{thin}{A positive integer specifying the period for saving samples. 
    The default is 1, which is usually the recommended value. Unless your
    posterior distribution takes up too much memory we do \emph{not} recommend 
    thinning as it throws away information. The tradition of thinning when 
    running MCMC stems primarily from the use of samplers that require 
    a large number of iterations to achieve the desired effective sample size. 
    Because of the efficiency (effective samples per second) of 
    Hamiltonian Monte Carlo, rarely should this be necessary when using Stan.}
    
  \item{init}{
    Specification of initial values for all or some parameters.
    Can be the digit \code{0}, the strings \code{"0"} or \code{"random"}, 
    a function that returns a named list, or a list of named lists:
    \describe{
    \item{\code{init="random"} (default):}{Let Stan generate random initial values
      for all parameters. The seed of the random number generator used by Stan 
      can be specified via the \code{seed} argument. If the seed for Stan is 
      fixed, the same initial values are used. The default is to randomly 
      generate initial values between \code{-2} and \code{2} \emph{on the 
      unconstrained support}. The optional additional parameter \code{init_r} can 
      be set to some value other than \code{2} to change the range of the randomly 
      generated inits.}
    \item{\code{init="0", init=0}:}{Initialize all parameters to zero on 
      the unconstrained support.}  
    \item{inits via list:}{Set inital values by providing a list 
      equal in length to the number of chains. The elements of this list should 
      themselves be named lists, where each of these named lists has the name 
      of a parameter and is used to specify the initial values for that parameter 
      for the corresponding chain.}
    \item{inits via function:}{Set initial values by providing a function that 
      returns a list for specifying the initial values of parameters for a chain. 
      The function can take an optional parameter \code{chain_id} through which the 
      \code{chain_id} (if specified) or the integers from 1 to \code{chains} will be 
      supplied to the function for generating initial values.
      See the \strong{Examples} section below for examples of defining 
      such functions and using a list of lists for specifying initial values.}
    }
    When specifying initial values via a \code{list} or \code{function}, any
    parameters for which values are not specified will receive initial values
    generated as described in the \code{init="random"} description above.
  }
  
  \item{seed}{The seed for random number generation. The default is generated 
    from 1 to the maximum integer supported by \R on the machine. Even if 
    multiple chains are used, only one seed is needed, with other chains having 
    seeds derived from that of the first chain to avoid dependent samples.
    When a seed is specified by a number, \code{as.integer} will be applied to it. 
    If \code{as.integer} produces \code{NA}, the seed is generated randomly. 
    The seed can also be specified as a character string of digits, such as
    \code{"12345"}, which is converted to integer.
    
    Using \R's \code{set.seed} function to set the seed for Stan will not work.
    }
    
  \item{algorithm}{
    One of the sampling algorithms that are implemented in Stan. 
    The default and preferred algorithm is \code{"NUTS"}, which is 
    the No-U-Turn sampler variant of Hamiltonian Monte Carlo
    (Hoffman and Gelman 2011, Betancourt 2017). Currently the other options 
    are \code{"HMC"} (Hamiltonian Monte Carlo), and \code{"Fixed_param"}. 
    When \code{"Fixed_param"} is used no MCMC sampling is performed 
    (e.g., for simulating with in the generated quantities block).
  }

  \item{sample_file}{An optional character string providing the name of a file.
    If specified the draws for \emph{all} parameters and other saved quantities
    will be written to the file. If not provided, files are not created. 
    When the folder specified is not writable, \code{tempdir()} is used. 
    When there are multiple chains, an underscore and chain number are appended
    to the file name.}

  \item{diagnostic_file}{An optional character string providing the name of a file.
    If specified the diagnostics data for \emph{all} parameters will be written
    to the file. If not provided, files are not created. When the folder specified 
    is not writable, \code{tempdir()} is used. When there are multiple chains, 
    an underscore and chain number are appended to the file name.}

  \item{save_dso}{Logical, with default \code{TRUE}, indicating whether the 
    dynamic shared object (DSO) compiled from the C++ code for the model 
    will be saved or not. If \code{TRUE}, we can draw samples from
    the same model in another \R session using the saved DSO (i.e., 
    without compiling the C++ code again).
    This parameter only takes effect if \code{fit} is not used; with  
    \code{fit} defined, the DSO from the previous run is used.
    When \code{save_dso=TRUE}, the fitted object can be loaded from 
    what is saved previously and used for sampling, if the compiling is
    done on the same platform, that is, same operating system and same 
    architecture (32bits or 64bits).}

  \item{verbose}{\code{TRUE} or \code{FALSE}: flag indicating whether 
    to print intermediate output from Stan on the console, which might
    be helpful for model debugging.}
 
  \item{control}{A named \code{list} of parameters to control the sampler's
    behavior. It defaults to \code{NULL} so all the default values are used. 
    First, the following are adaptation parameters for sampling algorithms.
    These are parameters used in Stan with similar names here.
    \itemize{
      \item \code{adapt_engaged} (\code{logical})
      \item \code{adapt_gamma} (\code{double}, positive, defaults to 0.05)
      \item \code{adapt_delta} (\code{double}, between 0 and 1, defaults to 0.8)
      \item \code{adapt_kappa} (\code{double}, positive, defaults to 0.75)
      \item \code{adapt_t0} (\code{double}, positive, defaults to 10)
      \item \code{adapt_init_buffer} (\code{integer}, positive, defaults to 75)
      \item \code{adapt_term_buffer} (\code{integer}, positive, defaults to 50)
      \item \code{adapt_window} (\code{integer}, positive, defaults to 25)
    }

    In addition, algorithm HMC (called 'static HMC' in Stan) and NUTS share the
    following parameters:
    \itemize{
      \item \code{stepsize} (\code{double}, positive)
      \item \code{stepsize_jitter} (\code{double}, [0,1])
      \item \code{metric} (\code{string}, one of "unit_e", "diag_e", "dense_e")
    }
    For algorithm NUTS, we can also set:
    \itemize{
      \item \code{max_treedepth} (\code{integer}, positive)
    }
    For algorithm HMC, we can also set:
    \itemize{
      \item \code{int_time} (\code{double}, positive)
    }
    
    For \code{test_grad} mode, the following parameters can be set:
    \itemize{
      \item \code{epsilon} (\code{double}, defaults to 1e-6)
      \item \code{error} (\code{double}, defaults to 1e-6)
    }
  }
  
  \item{open_progress}{Logical scalar that only takes effect if 
    \code{cores > 1} but is recommended to be \code{TRUE} in interactive
    use so that the progress of the chains will be redirected to a file
    that is automatically opened for inspection. For very short runs, the
    user might prefer \code{FALSE}.}
  
  \item{\dots}{Other optional parameters: 
    \itemize{
      \item \code{chain_id} (\code{integer})
      \item \code{init_r} (\code{double}, positive)
      \item \code{test_grad} (\code{logical})
      \item \code{append_samples} (\code{logical})
      \item \code{refresh}(\code{integer})
      \item \code{save_warmup}(\code{logical})
      \item deprecated: \code{enable_random_init}(\code{logical}) 
    }

    \code{chain_id} can be a vector to specify the chain_id for all 
    chains or an integer. For the former case, they should be unique. 
    For the latter, the sequence of integers starting from the given
    \code{chain_id} are used for all chains. 

    \code{init_r} is used only for generating random initial values,
    specifically when \code{init="random"} or not all parameters
    are initialized in the user-supplied list or function. If specified, 
    the initial values are simulated uniformly from interval 
    [-\code{init_r}, \code{init_r}] rather than using the default interval 
    (see the manual of (cmd)Stan).
    
    \code{test_grad} (\code{logical}).  
    If \code{test_grad=TRUE}, Stan will not do any sampling. Instead,
    the gradient calculation is tested and printed out and the fitted
    \code{stanfit} object is in test gradient mode.  By default, it is
    \code{FALSE}.
    
    \code{append_samples} (\code{logical}). 
    Only relevant if \code{sample_file} is specified \emph{and} is an 
    existing file. In that case, setting \code{append_samples=TRUE} will append 
    the samples to the existing file rather than overwriting the contents of 
    the file.

    \code{refresh} (\code{integer}) can be used to
    control how often the progress of the sampling is reported (i.e. 
    show the progress every \code{refresh} iterations). 
    By default, \code{refresh = max(iter/10, 1)}. 
    The progress indicator is turned off if \code{refresh <= 0}.

    Deprecated: \code{enable_random_init} (\code{logical}) being \code{TRUE}
    enables specifying initial values randomly when the initial
    values are not fully specified from the user. 
    
    \code{save_warmup} (\code{logical}) indicates whether to 
    save draws during the warmup phase and defaults to \code{TRUE}. Some
    memory related problems can be avoided by setting it to \code{FALSE},
    but some diagnostics are more limited if the warmup draws are not
    stored.
  }

  \item{boost_lib}{The path for an alternative version of the Boost C++ 
    to use instead of the one in the \pkg{BH} package.}

  \item{eigen_lib}{The path for an alternative version of the Eigen C++ 
     library to the one in \pkg{RcppEigen}.} 
}
\details{
  The \code{stan} function does all of the work of fitting a Stan model and 
  returning the results as an instance of \code{stanfit}. The steps are 
  roughly as follows:
  
  \enumerate{
  \item Translate the Stan model to C++ code.  (\code{\link{stanc}})
  \item Compile the C++ code into a binary shared object, which
  is loaded into the current \R session (an object
  of S4 class \code{stanmodel} is created).   (\code{\link{stan_model}})
  \item Draw samples and wrap them in an object of S4 class \code{stanfit}. (\code{\link{sampling}}) 
  }
  
  The returned object can be used with methods such as \code{print}, 
  \code{summary}, and \code{plot} to inspect and retrieve the results of 
  the fitted model. 

  \code{stan} can also be used to sample again from a fitted model under
  different settings (e.g., different \code{iter}, \code{data}, etc.) by 
  using the \code{fit} argument to specify an existing \code{stanfit} object. 
  In this case, the compiled C++ code for the model is reused. 
}
\value{
  An object of S4 class \code{\link{stanfit}}. However, if \code{cores > 1}
  and there is an error for any of the chains, then the error(s) are printed. If 
  all chains have errors and an error occurs before or during sampling, the returned
  object does not contain samples. But the compiled binary object for the 
  model is still included, so we can reuse the returned object for another
  sampling. 
}

\section{Passing data to Stan}{
  The data passed to \code{stan} are preprocessed before being passed to Stan.
  If \code{data} is not a character vector, the data block of the Stan program
  is parsed and \R objects of the same name are searched starting from the
  calling environment. Then, if \code{data} is list-like but not a \code{data.frame}
  the elements of \code{data} take precedence. This behavior is similar to how
  a \code{formula} is evaluated by the \code{\link{lm}} function when \code{data} is
  supplied. In general, each \R object being passed to Stan should be either a numeric
  vector (including the special case of a 'scalar') or a numeric array (matrix). 
  The first exception is that an element can be a logical vector: \code{TRUE}'s 
  are converted to 1 and \code{FALSE}'s to 0.
  An element can also be a data frame or a specially structured list (see
  details below), both of which will be converted into arrays in the
  preprocessing.  Using a specially structured list is not
  encouraged though it might be convenient sometimes; and when in doubt, just
  use arrays.  

  This preprocessing for each element mainly includes the following:
  \enumerate{
    \item Change the data of type from \code{double} 
      to \code{integer} if no accuracy is lost. The main
      reason is that by default, \R uses \code{double}
      as data type such as in \code{a <- 3}. But Stan 
      will not read data of type \code{int} from \code{real} 
      and it reads data from \code{int} if the data 
      type is declared as \code{real}. 
   
   \item Check if there is \code{NA} in the data. 
     Unlike BUGS, Stan does not allow missing data. Any \code{NA} values
     in supplied data will cause the function to stop and report an error.
  
   \item Check data types.  Stan allows only numeric data, that is,
     doubles, integers, and arrays of these.  Data of other types (for
     example, characters and factors) are not passed to Stan. 

   \item Check whether there are objects in the data list with duplicated names.
     Duplicated names, if found, will cause the function to stop and
     report an error. 

   \item Check whether the names of objects in the data list are legal 
     Stan names. If illegal names are found, it will stop and 
     report an error. See (Cmd)Stan's manual for the rules of variable names. 

   \item When an element is of type \code{data.frame}, it will be converted to
     \code{matrix} by function \code{data.matrix}. 

   \item When an element is of type \code{list}, it is supposed to make it
     easier to pass data for those declared in Stan code such as 
     \code{"vector[J] y1[I]"} and \code{"matrix[J,K] y2[I]"}. Using the latter
     as an example, we can use a list for \code{y2} if the list has "I" elements,
     each of which is an array (matrix) of dimension "J*K". However, it is 
     not possible to pass a list for data declared such as 
     \code{"vector[K] y3[I,J]"}; the only way for it is to use an array with
     dimension "I*J*K". In addition, technically a \code{data.frame} in R is 
     also a list, but it should not be used for the purpose here since a
     \code{data.frame} will be converted to a matrix as described above.  
  } 

  Stan treats a vector of length 1 in R as a scalar.  So technically 
  if, for example, \code{"real y[1];"} is defined in the data block, an array 
  such as \code{"y = array(1.0, dim = 1)"} in R should be used. This
  is also the case for specifying initial values since the same 
  underlying approach for reading data from R in Stan is used, in which
  vector of length 1 is treated as a scalar.

  In general, the higher the optimization level is set, the faster the 
  generated binary code for the model runs, which can be set in a Makevars file.
  However, the binary code generated for the model runs fast by using a higher
  optimization level at the cost of longer times to compile the C++ code. 
}

\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{http://mc-stan.org}. 

  The Stan Development Team
  \emph{CmdStan Interface User's Guide}.
  \url{http://mc-stan.org}.
}
\seealso{
  \itemize{
  \item The package vignettes for an example of fitting a model and accessing
  the contents of \code{stanfit} objects (\url{http://mc-stan.org/rstan/articles/}).
  
  \item \code{\link{stanc}} for translating model code in Stan modeling language to C++, 
  \code{\link{sampling}} for sampling, and \code{\linkS4class{stanfit}} for the
  fitted results.

  \item \code{\link{as.array.stanfit}} and \code{\link{extract}} for extracting
  samples from \code{stanfit} objects.
  }
}
\examples{\dontrun{
#### example 1 
library(rstan)
scode <- "
parameters {
  real y[2]; 
} 
model {
  y[1] ~ normal(0, 1);
  y[2] ~ double_exponential(0, 2);
} 
"
fit1 <- stan(model_code = scode, iter = 10, verbose = FALSE) 
print(fit1)
fit2 <- stan(fit = fit1, iter = 10000, verbose = FALSE) 

## using as.array on the stanfit object to get samples 
a2 <- as.array(fit2)

## extract samples as a list of arrays
e2 <- extract(fit2, permuted = FALSE)

#### example 2
#### the result of this package is included in the package 

excode <- '
  transformed data {
    real y[20];
    y[1] = 0.5796;  y[2] = 0.2276;   y[3]  = -0.2959; 
    y[4] = -0.3742; y[5] = 0.3885;   y[6]  = -2.1585;
    y[7] = 0.7111;  y[8] = 1.4424;   y[9]  = 2.5430; 
    y[10] = 0.3746; y[11] = 0.4773;  y[12] = 0.1803; 
    y[13] = 0.5215; y[14] = -1.6044; y[15] = -0.6703; 
    y[16] = 0.9459; y[17] = -0.382;  y[18] = 0.7619;
    y[19] = 0.1006; y[20] = -1.7461;
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

exfit <- stan(model_code = excode, save_dso = FALSE, iter = 500)
print(exfit)
plot(exfit)
}
\dontrun{
## examples of specify argument `init` for function stan

## define a function to generate initial values that can
## be fed to function stan's argument `init`
# function form 1 without arguments 
initf1 <- function() {
  list(mu = 1, sigma = 4, z = array(rnorm(6), dim = c(3,2)), alpha = 1)
} 
# function form 2 with an argument named `chain_id`
initf2 <- function(chain_id = 1) {
  # cat("chain_id =", chain_id, "\n")
  list(mu = 1, sigma = 4, z = array(rnorm(6), dim = c(3,2)), alpha = chain_id)
} 

# generate a list of lists to specify initial values
n_chains <- 4
init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))
 
exfit0 <- stan(model_code = excode, init = initf1) 
stan(fit = exfit0, init = initf2) 
stan(fit = exfit0, init = init_ll, chains = n_chains)
}} 
\keyword{rstan} 
