# This file is part of RStan
# Copyright (C) 2015 Jiqiang Guo and Benjamin Goodrich
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
#' header files of the Stan C++ library.
#' 
#' Starting from \pkg{rstan} 3.0, the syntax of \code{\link{ReferenceClasses}}
#' is heavily utilized by the \pkg{rstan} package. This syntax is unfamiliar
#' to many \R users but is similar to that in Python. In particular, various
#' \code{\link{ReferenceClasses}} are defined with fields and methods. These
#' methods are \code{\link{function}}s that may utilize the fields in addition
#' to any arguments of the methods. Perhaps the biggest syntactic difference is
#' that these methods are called using the \code{$} operator, such as
#' \code{object$plot()} rather than \code{plot(object)}. To list the methods
#' that are defined for \code{\link{ReferenceClasses}}, call the 
#' \code{$methods()} method with no arguments or simply postfix \code{$} to
#' the object name and press the Tab key. To obtain help for any defined
#' method, use the \code{$help(topic)} method, where \code{topic} is the name
#' of the method you need help with.
#' 
#' The steps for defining, estimating, and diagnosing a Stan model are as
#' follows:
#' 
#' @section Step 0 --- write a Stan program:
#' Preferably save a file with a .stan extension to your hard disk that
#' constitutes a valid Stan program or alternatively you can write a valid
#' Stan program as a \R \code{\link{character}} vector, in which case it will
#' be copied into a .stan file in the \code{\link{tempdir}} whose prefix is its
#' \code{\link[tools]{md5sum}}. A simple illustration of a Stan program is
#' given in the examples section below. More details about the Stan language
#' are at \url{http://mc-stan.org/manual.html}, but in brief a Stan program has
#' one or more of the following blocks contained within curly braces:
#' \itemize{
#'   \item functions: an optional section for users to define functions that
#'     can subsequently be used in that Stan program only or can be exported
#'     to the \R \code{\link{.GlobalEnv}}
#'   \item data: an all-but-required section for users to declare known 
#'     quantities, including the outcome to be modeled, its predictors,
#'     hyperparameters, dimensions, etc.
#'   \item transformed data: an optional section for users to define other
#'     known quantities that are deterministically related to the objects
#'     declared in the \code{data} block and hence is only executed once
#'   \item parameters: an all-but-required section to declare the unknowns
#'     that are being estimated
#'   \item transformed parameters: an optional section for users to define
#'     unknown quantities that are conditionally deterministic given the
#'     objects in \code{parameters}, \code{transformed data}, and \code{data}.
#'     The \code{transformed parameters} block is estimated every time the
#'     objective function is evaluated, the constraints on the objects therein
#'     are checked, and the objects are saved in the output
#'   \item model: an all-but-required section for users to build up an
#'     expression for the objective function, which may consist of both priors
#'     and a likelihood to form a posterior. Local objects can be declared at
#'     the top of the \code{model} block and subsequently defined but their
#'     constraints are not validated and they are not saved in the output.
#'   \item generated quantities: an optional section where users can define
#'     additional objects that depend on anything in the previous sections
#'     \emph{except} the \code{model} section and may involve random-number
#'     generation. The \code{generated quantities} section is executed once
#'     per iteration and is useful for producing posterior predictive 
#'     distributions.
#' }
#' Stan supports the several types of objects. The primitive scalar types are
#' \itemize{
#'   \item \code{real} which can be unknown or known
#'   \item \code{int} although unknown quanties cannot be integers  
#' }
#' The unconstrained linear algebra types contain real numbers and are
#' \itemize{
#'   \item \code{vector[K]} a (column) vector of length K
#'   \item \code{row_vector[K]} a row vector of length K
#'   \item \code{matrix[M,N]} a matrix with M rows and N columns even if M or
#'     N is 1
#' }
#' Any primitive or unconstrained linear algebra type can optionally be
#' inequality-constrained by appending \code{<lower=l,upper=u>} to the
#' last letter of the object type (i.e. before the [ in the case of a linear
#' algebra type). The scalar values \code{l} and \code{u} may be expressions
#' of previously declared objects, may be infinite, and either may be omitted,
#' in which case they are taken to be infinite.
#' The constrained linear algebra types contain real numbers, cannot be
#' declared within the \code{model} block or within the body of a user-defined
#' function, cannot be further constrained by inequalities, and are
#' \itemize{
#'   \item \code{simplex[K]} a non-negative (column) vector of length K that
#'     sums to 1
#'   \item \code{unit_vector[K]} a (column) vector of length K whose squares
#'     sum to 1
#'   \item \code{cov_matrix[K]} a square, positive-definite matrix of order K
#'   \item \code{corr_matrix[K]} a square, positive-definite matrix of order K
#'     whose diagonal elements are all 1
#'   \item \code{cholesky_factor_cov[M,N]} a matrix with M rows and N <= M 
#'     columns whose diagonal elements are positive and whose above-diagonal
#'     elements are zero
#'   \item \code{cholesky_factor_corr[K]} a square matrix of order K whose
#'     diagonal elements are positive, whose above-diagonal elements are zero,
#'     and whose rows are unit vectors 
#' }
#' Finally, the user can declare an array of \emph{any} of the above objects by
#' postfixing its dimensions within square brackets after the object's name.
#' For example, \code{int<lower=0> y[N];} would declare an array of length 
#' \code{N} named \code{y} that consists of non-negative integers. An array in
#' Stan is similar to an \R list but whose elements are homogenous. Also,
#' arrays can be multidimensional with the dimensions separated by commas, 
#' which would be like a \R list of lists.
#' 
#' @section Step 1 --- Create an object of \code{\link{StanProgram-class}}:
#' Either pass a file path to the \code{file} argument of the 
#' \code{StanProgram} constructor or pass the model in \code{\link{character}}
#' form to the \code{argument} of the \code{StanProgram} constructor. For
#' example, \code{program <- StanProgram(file = MyStanFile.stan)} would locate
#' a file called MyStanFile.stan in the working directory on the hard disk. An
#' illustration of the \code{model} argument is given in the examples section
#' below. Either approach will cause the Stan syntax to be parsed to C++ and
#' if successful compiled, which may take a considerable amount of time. There 
#' are three things that users would often want to do with such an object:
#' \itemize{
#'   \item save the object to the hard disk via the \code{$save()} method,
#'     which eliminates the need to compile the same Stan program again (even
#'     in a different \R session if the \code{file} argument was specified)
#'     and is invoked automatically if \code{rstan_options(auto_write = TRUE)}
#'     is first executed
#'   \item expose any user-defined Stan functions via the \code{$expose()} 
#'     method, which exposes them to the \R \code{\link{.GlobalEnv}} in order
#'     to verify that they are working correctly or to use as part of a
#'     simulation
#'   \item instantiate an object of \code{\link{StanProgramWithData-class}} via
#'     the \code{$instantiate()} method which associates the \code{StanProgram}
#'     with the \R objects that correspond to the objects declared in the
#'     \code{data} block of the Stan program
#' }
#' 
#' @section Step 2 --- Create a \code{\link{StanProgramWithData-class}} object:
#' In particular, it is eventually necessary to instantiate an object of
#' \code{\link{StanProgramWithData-class}} in order to estimate the parameters.
#' For example, a user could execute \code{dprogram <- program$instantiate()}.
#' The \code{$instantiate(data)} method takes an optional \code{data} argument,
#' which can be a named list of \R objects, an environment containing \R 
#' objects, or a \code{\link{data.frame}} that is passed to the 
#' \code{\link{data.matrix}} function to create objects from the resulting
#' columns. However, if the \code{data} argument is omitted, then the calling
#' \code{\link{environment}} is searched for \R objects whose symbol names are
#' the same as the names of the objects declared in the \code{data} block of
#' the Stan program, which is usually the most convenient way to proceed.
#' 
#' There are two things that users would often want to do with such an object:
#' \itemize{
#'   \item \code{$sample()} which performs Markov Chain Monte Carlo (MCMC)
#'     using the default sampling algorithm and outputs an object of
#'     \code{\link{StanFitMCMC-class}}
#'   \item \code{$optimize()} which finds the posterior mode using the default
#'     optimization algorithm and outputs an object of
#'     \code{\link{StanFitOptimization-class}}
#' }
#' However, it may be essential or desireable to use a non-default sampling
#' or optimization algorithm or to specify non-default arguments. In that case,
#' it is necessary to call one of the following methods:
#' \itemize{
#'   \item sampling algorithms
#'   \itemize{
#'     \item ehmc which is the default sampling algorithm but accepts
#'       additional arguments that may be given non-default values
#'     \item MORE
#'   }
#'   \item optimizaiton algorithms
#'   \itemize{
#'     \item lfbgs which is the default optimization algorithm but accepts
#'       additional arguments that may be given non-default values
#'     \item bfgs
#'     \item newton
#'   }
#' }
#' The sampling algorithms take an overlapping set of optional arguments, but
#' the most common arguments that a user might want to specify non-default
#' values for are
#' \itemize{
#'   \item \code{chains} which indicates how many Markov chains to execute
#'   \item \code{cores} which indicates how many of the local computer's 
#'     cores to utilize and should set equal to the number of chains if
#'     RAM permits
#'   \item \code{delta} which indicates the target acceptance rate and
#'     should be set to a value higher than the default of 0.8 if the
#'     default results in any divergent transitions after the warmup
#'     phase (see Step X below)
#'   \item MORE
#' }
#' The optimization algorithms take an overlapping set of optional arguments,
#' but the most common arguments that a user might want to specify non-default
#' values for are
#' \itemize{
#'   \item \code{rel_tol}
#'   \item \code{abs_tol}
#'   \item MORE
#' }
#' In addition, there are many more utility methods that some users may wish to
#' call that are documented in \code{\link{StanProgramWithData-methods}}.
#' 
#' @section Step 3 --- Estimate the parameters:
#' Call some estimation method and assign the result to an object of
#' \code{\link{StanFitMCMC-class}} in the case of a sampling method or
#' \code{\link{StanFitOptimization-class}} in the case of an optimization
#' method. For example, execute \code{posterior <- dprogram$sample()}.
#' 
#' @section Step 4 --- Diagnose any problems:
#' To be written
#' 
#' @examples
#' # create data in R
#' J <- 8                                     # number of schools
#' y <- c(28,  8, -3,  7, -1,  1, 18, 12)     # estimated treatment effects
#' sigma <- c(15, 10, 16, 11,  9, 11, 10, 18) # estimated standard deviations
#' 
#' # Step 0 --- write a Stan program 
#' mc <- '
#' data { /* objects conditioned on in Bayes' theorem */
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
#'   vector[J] theta;
#'   theta <- mu + tau * eta;
#' }
#' model { /* build up expression for the posterior distribution */
#'   # implicit: improper priors on mu and tau
#'   eta ~ normal(0, 1);       // normal prior on eta
#'   # implies: theta ~ normal(mu, tau)
#'   y ~ normal(theta, sigma); // normal likelihood for y | theta
#' }
#' '
#' 
#' \donttest{
#' # Step 1 --- Create an object of StanProgram-class
#' program <- StanProgram(code = mc) # preferable to specify a file
#'                                           
#' # Step 2 --- Create a StanProgramWithData-class object
#' dprogram <- program$instantiate()
#' 
#' # Step 3 --- Estimate the parameters
#' estimates <- dprogram$optimize()        # maximum a posteriori estimator
#' estimates <- dprogram$ehmc(delta = 0.9) # MCMC from posterior distribution
#' 
#' # Step 4 --- Diagnose any problems
#' estimates$summary()
#' }
#' @docType package
#' @name rstan3
NULL