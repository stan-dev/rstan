\name{sampling-methods}
\docType{methods}
\alias{sampling} 
\alias{sampling-methods}
\alias{sampling,stanmodel-method}

\title{Draw samples from a Stan model}
\description{
  Draw samples from the model defined by class \code{stanmodel}. 
  This method is a generic function of the S4 class \code{stanmodel}. 
}

\usage{
  %% sampling(object, \dots)  
  \S4method{sampling}{stanmodel}(object, data = list(), pars = NA, 
    chains = 4, iter = 2000, warmup = floor(iter/2), thin = 1,
    seed = sample.int(.Machine$integer.max, 1), 
    init = 'random', check_data = TRUE, 
    sample_file = NULL, diagnostic_file = NULL, verbose = FALSE, 
    algorithm = c("NUTS", "HMC", "Fixed_param"),
    control = NULL, include = TRUE, 
    cores = getOption("mc.cores", 1L),
    open_progress = interactive() && !isatty(stdout()) &&
                    !identical(Sys.getenv("RSTUDIO"), "1"),
    show_messages = TRUE, \dots)
} 

\section{Methods}{
  \describe{
    \item{\code{sampling}}{\code{signature(object = "stanmodel")}}{
      %%  ~~describe this method here~~
      Call a sampler (NUTS, HMC, or Fixed_param depending on parameters) 
      to draw samples from the model defined by S4 class \code{stanmodel} 
      given the data, initial values, etc. 
    }
  }
}

\arguments{
  \item{object}{An object of class \code{\linkS4class{stanmodel}}.}

  \item{data}{A named \code{list} or \code{environment}
    providing the data for the model or a character vector 
    for all the names of objects used as data. 
    See the Note section in \code{\link{stan}}.} 

  \item{pars}{A vector of character strings specifying parameters of interest. 
    The default is \code{NA} indicating all parameters in the model. 
    If \code{include = TRUE}, only samples for parameters named in \code{pars} 
    are stored in the fitted results. Conversely, if \code{include = FALSE}, 
    samples for all parameters \emph{except} those named in \code{pars} are 
    stored in the fitted results.}

  \item{chains}{A positive integer specifying the number of Markov chains. 
    The default is 4.} 
 
  \item{iter}{A positive integer specifying the number of iterations for each 
    chain (including warmup). The default is 2000.}

  \item{warmup}{A positive integer specifying the number of warmup (aka burnin)
    iterations per chain. If step-size adaptation is on (which it is by default), 
    this also controls the number of iterations for which adaptation is run (and
    hence these warmup samples should not be used for inference). The number of 
    warmup iterations should not be larger than \code{iter} and the default is 
    \code{iter/2}.}

  \item{thin}{A positive integer specifying the period for saving samples. 
    The default is 1, which is usually the recommended value.}

  \item{seed}{The seed for random number generation. The default is generated 
    from 1 to the maximum integer supported by \R on the machine. Even if 
    multiple chains are used, only one seed is needed, with other chains having 
    seeds derived from that of the first chain to avoid dependent samples.
    When a seed is specified by a number, \code{as.integer} will be applied to it. 
    If \code{as.integer} produces \code{NA}, the seed is generated randomly. 
    The seed can also be specified as a character string of digits, such as
    \code{"12345"}, which is converted to integer.}

  \item{init}{Initial values specification. See the detailed documentation for 
    the init argument in \code{\link{stan}}.}

  \item{check_data}{Logical, defaulting to \code{TRUE}. If \code{TRUE} 
    the data will be preprocessed; otherwise not.
    See the Note section in \code{\link{stan}}.}

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
   
  \item{verbose}{\code{TRUE} or \code{FALSE}: flag indicating whether 
    to print intermediate output from Stan on the console, which might
    be helpful for model debugging.}

  \item{algorithm}{One of algorithms that are implemented in Stan such 
    as the No-U-Turn sampler (NUTS, Hoffman and Gelman 2011), static HMC,
    or Fixed_param.} 

  \item{control}{A named \code{list} of parameters to control the sampler's
    behavior. See the details in the documentation for the \code{control} argument
    in \code{\link{stan}}.}

  \item{include}{Logical scalar defaulting to \code{TRUE} indicating
    whether to include or exclude the parameters given by the 
    \code{pars} argument. If \code{FALSE}, only entire multidimensional
    parameters can be excluded, rather than particular elements of them.}
  
  \item{cores}{Number of cores to use when executing the chains in parallel,
    which defaults to 1 but we recommend setting the \code{mc.cores} option 
    to be as many processors as the hardware and RAM allow (up to the 
    number of chains).}
    
  \item{open_progress}{Logical scalar that only takes effect if 
    \code{cores > 1} but is recommended to be \code{TRUE} in interactive
    use so that the progress of the chains will be redirected to a file
    that is automatically opened for inspection. For very short runs, the
    user might prefer \code{FALSE}.}
    
  \item{show_messages}{Either a logical scalar (defaulting to \code{TRUE})
    indicating whether to print the summary of Informational Messages to
    the screen after a chain is finished or a character string naming a path
    where the summary is stored. Setting to \code{FALSE} is not recommended
    unless you are very sure that the model is correct up to numerical 
    error.}

  \item{\dots}{Additional arguments can be \code{chain_id}, \code{init_r}, 
    \code{test_grad}, \code{append_samples}, \code{refresh},
    \code{enable_random_init}. See the documentation in \code{\link{stan}}.} 
}

\value{
   An object of S4 class \code{stanfit} representing
   the fitted results. Slot \code{mode} for this object
   indicates if the sampling is done or not. 
} 

\seealso{
  \code{\linkS4class{stanmodel}}, \code{\linkS4class{stanfit}}, \code{\link{stan}}
}

\examples{\dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
f <- sampling(m, iter = 100)
}}
