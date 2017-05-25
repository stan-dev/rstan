# This file is part of RStan
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
#
# RStan is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# RStan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Estimate and diagnose models using Stan
#' 
#' This R package provides an interface to the Stan C++ library or more
#' specifically links to the \pkg{StanHeaders} package, which provides the
#' header files of the Stan C++ libraries.
#' 
#' Starting from \pkg{rstan} 3.0, the syntax of \code{\link[methods]{ReferenceClasses}}
#' is heavily utilized by the \pkg{rstan} package. This syntax may be unfamiliar
#' to many \R users but is similar to that in Python. In particular, various
#' \code{\link{ReferenceClasses}} are defined with fields and methods. These
#' methods are \code{\link{function}}s that may utilize the fields in addition
#' to any arguments of the methods. Perhaps the biggest syntactic difference is
#' that these methods are called using the \code{$} operator, such as
#' \code{object$plot()} rather than \code{plot(object)}. To list the methods
#' that are defined for \code{\link[methods]{ReferenceClasses}}, postfix \code{$}
#' to the object name and press the Tab key. To obtain help for any defined
#' method, use the \code{$help(topic)} method, where \code{topic} is the name
#' of the method you need help with.
#' 
#' The basic steps for defining, estimating, and diagnosing a Stan model are as
#' follows.
#' 
#' @section Step 0 --- Write a Stan program:
#' Preferably save a file with a .stan extension to your hard disk that
#' constitutes a valid Stan program or you can write a valid
#' Stan program as a \R \code{\link{character}} vector. A simple illustration 
#' of a Stan program is given in the examples section below. More details about 
#' the Stan language are at \url{http://mc-stan.org/documentation}, but in brief 
#' a Stan program has one or more of the following blocks contained within curly 
#' braces:
#' \itemize{
#'   \item functions: an optional section for users to define functions that
#'     can subsequently be used in that Stan program only or can be exported
#'     to the \R by calling \code{link[rstan]{expose_stan_functions}}
#'   \item data: an all-but-required section for users to declare known 
#'     quantities, including the outcome to be modeled, its predictors,
#'     hyperparameters, dimensions, etc. These are passed from R.
#'   \item transformed data: an optional section for users to define other
#'     known quantities that are related to the objects declared in the \code{data} 
#'     block and hence is only executed once. It is also possible to do random
#'     number generation in this block.
#'   \item parameters: an all-but-required section to declare the unknowns
#'     that are being estimated and hence are saved in the output
#'   \item transformed parameters: an optional section for users to define
#'     unknown quantities that are conditionally deterministic given the
#'     objects in \code{parameters}, \code{transformed data}, and \code{data}.
#'     The \code{transformed parameters} block is executed every time the
#'     objective function is evaluated, the constraints on the objects therein
#'     are checked, and the objects are saved in the output
#'   \item model: an all-but-required section for users to build up an
#'     expression for the objective function, which may consist of both priors
#'     and a likelihood to form a posterior kernel. Local objects can be declared
#'     at the top of the \code{model} block but their are not saved in the output.
#'   \item generated quantities: an optional section where users can define
#'     additional objects that depend on anything in the previous sections
#'     \emph{except} the \code{model} section and may involve random-number
#'     generation. The \code{generated quantities} section is executed once
#'     per iteration and is useful for producing posterior predictive 
#'     distributions.
#' }
#' The Stan language supports the several types of objects. These same types are now
#' available in \R via S4 classes, as documented by the \code{\link{StanType-class}}.
#' 
#' @section Step 1 --- Compile your Stan program:
#' Either pass a file path to the \code{file} argument of the 
#' \code{\link{stan_config}} constructor or pass specify the \code{mode_code}
#' argument. For example, \code{config <- stan_config(file = "MyStanFile.stan")} 
#' would locate a file called MyStanFile.stan in the working directory on the 
#' hard disk. If the \code{file} argument is unspecified, a dialogue box will
#' open (in \code{\link[base]{interactive}} mode) that will prompt you to a file
#' on your hard disk. An illustration of the \code{model_code} argument is given  
#' in the examples section below. Either approach will cause the Stan syntax to
#' be parsed to C++ and (if successful) compiled, which may take a considerable 
#' amount of time. This will create an object of \code{\link{StanConfig-class}}.
#' 
#' @section Step 2 --- Prepare to estimate the parameters of your model:
#' First, specify what is being passed from \R to the \code{data} block of the
#' Stan program using the \code{$data} method for an object or 
#' \code{\link{StanConfig-class}}.These \R objects can either be specified by 
#' name or the \code{$data} methodcan be called without any arguments to 
#' \code{\link[base]{dynGet}} them from the inherited 
#' \code{\link[base]{environment}}(s) that the object of
#' \code{\link{StanProgram-class}} was created in.
#' 
#' Second, once the objects in the \code{data} block of the Stan program have
#' been fully specified, call the \code{$build} method to pass them to the C++
#' constructor and instantiate a Stan model object, which in \R is represented
#' as a \pkg{Rcpp} \code{\link[Rcpp]{Module}} that is contained within an object
#' of \code{\link{StanModule-class}}.
#' 
#' @section Step 3 --- Estimate the unknowns your model:
#' There are a variety of methods for an object of \code{\link{StanModule-class}}
#' that can be invoked to estimate the unknowns of your model in various ways.
#' The preferred way is to draw from the posterior distribution using Markov 
#' Chain Monte Carlo (MCMC) sampling, but it is also possible to use Automatic 
#' Differentiation Variational Inference (ADVI) with either a meanfield or a 
#' fullrank assumption on the covariance matrix of the multivariate normal 
#' distribution that is closest to the posterior distribution in the 
#' unconstrained space. Finally, it is possible to find the mode of a
#' log-likelihood (if not priors are specified) or a posterior mode (if priors
#' are specified) using various deterministic optimization algorithms that also
#' make use of automatic differentiation to evaluate the gradient with respect
#' to the unknowns in the unconstrained parameter space. These algorithms can
#' be called with their default arguments by invoking one of the following
#' methods:
#' \itemize{
#'   \item \code{sample}: MCMC, specifically No U-Turn Sampling (NUTS) with
#'     adaptive integration time, a tuned diagonal mass matrix, and a 
#'     Euclidean metric
#'   \item \code{advi}: ADVI with a diagonal covariance matrix
#'   \item \code{maximize}: The LBFGS mode-finding algorithm
#' } 
#' Other algorithms can be invoked (possibly with non-default arguments) via 
#' the corresponding methods, which are documented with the 
#' \code{\link{StanModule-class}}.
#'  
#' @section Step 4 --- Diagnose any problems:
#' The best way to diagnose problems with MCMC is to call the 
#' \code{\link[shinyStan]{launch_shinystan}} function in the 
#' \pkg{shinyStan} package, but this requires user interactivity. Many
#' static plots can be produced interactively or by a script (via the 
#' \pkg{bayesplot} package) by invoking the methods for an object of 
#' \code{\link{StanMCMCFit-class}}. However, the most basic methods are 
#' \code{show} and, for more detail, \code{summary}.
#' 
#' @examples
#' # create data in R
#' J <- 8                                     # number of schools
#' y <- c(28,  8, -3,  7, -1,  1, 18, 12)     # estimated treatment effects
#' sigma <- c(15, 10, 16, 11,  9, 11, 10, 18) # estimated standard deviations
#' 
#' # Step 0 --- Write a Stan program 
#' mc <- '
#' data { /* objects conditioned on in Bayes Law */
#'   int<lower=0> J;           // syntax for an non-negative integer scalar 
#'   real y[J];                // syntax for an array of real numbers
#'   vector<lower=0>[J] sigma; // syntax for a non-negative vector of reals
#' }
#' parameters { /* variables whose posterior distribution is sought */
#'   real mu; 
#'   real<lower=0> tau;
#'   vector[J] eta;
#' }
#' transformed parameters { /* conditionally deterministic variables */
#'   vector[J] theta = mu + tau * eta;
#' }
#' model { /* build up expression for the posterior distribution */
#'   // implicit: improper priors on mu and tau
#'   target += normal_lpdf(eta | 0, 1);
#'   // normal prior on eta implies: theta ~ normal(mu, tau)
#'   target += normal_lpdf(y | theta, sigma); // normal likelihood
#' }
#' '
#' \donttest{
#' # Step 1 --- Compile your Stan program
#' # a path to file is preferable but a character works
#' config <- stan_config(model_code = mc, model_name = "8schools")
#'                                           
#' # Step 2 --- Prepare to estimate the parameters of your model
#' config$data()            # grabs J, y, and sigma from .GlobalEnv
#' config$data(J, y, sigma) # has the same effect as the previous line
#' module <- config$build(seed = 12345)
#' 
#' # Step 3 --- Estimate the unknowns ofyour model
#' module$help("hmc_nuts_diag_e_adapt") # this is what is called by $sample()
#' post <- module$sample()              # draws from posterior distribution
#' 
#' # Step 4 --- Diagnose any problems
#' post           # shows basic characteristics of the posterior distribution
#' post$summary() # shows more details about the posterior distribution
#' }
#' @docType package
#' @name rstan3
NULL
