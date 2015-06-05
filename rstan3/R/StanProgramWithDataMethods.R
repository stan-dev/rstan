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

StanProgramWithData$methods(help = help_from_instance)

StanProgramWithData$methods(log_prob = function(uparams) {
  "Value of the log-posterior at the unconstrained parameters"
  stop("FIXME: Implement")
})

StanProgramWithData$methods(grad = function(uparams) {
  "Gradient of the log-posterior at the unconstrained parameters"
  stop("FIXME: Implement")
})

StanProgramWithData$methods(log_prob_grad = function(uparams) {
  "Returns list with the value and gradient of the log-posterior at the unconstrained parameters"
  stop("FIXME: Implement")
})

StanProgramWithData$methods(hessian = function(uparams) {
  "Hessian of the log-posterior at the unconstrained parameters"
  stop("FIXME: Implement")
})

StanProgramWithData$methods(log_prob_grad_hessian = function(uparams) {
  "Returns list with the value, gradient, and Hessian of the log-posterior at the unconstrained parameters"
  stop("FIXME: Implement") # Does not exist in C++ yet
})

StanProgramWithData$methods(laplace_approx = function(uparams) {
  "Laplace approximation of the normalizing factor expanded at the unconstrained parameters"
  stop("FIXME: Implement") # Does not exist in C++ yet
})

StanProgramWithData$methods(constrain_params = function(uparams) {
  "Returns list of constrained parameters after transformation from the unconstrained parameters"
  stop("FIXME: Implement")
})

StanProgramWithData$methods(unconstrain_params = function(params) {
  "Returns vector of unconstrained parameters after transformation from the constrained parameters"
  stop("FIXME: Implement")
})

StanProgramWithData$methods(lbfgs = function() {
  "Find posterior mode with the LBFGS algorithm"
  stop("FIXME: Implement")
})

StanProgramWithData$methods(bfgs = function() {
  "Find posterior mode with the BFGS algorithm"
  stop("FIXME: Implement")
})

StanProgramWithData$methods(newton = function() {
  "Find posterior mode with the Newton algorithm"
  stop("FIXME: Implement")
})

StanProgramWithData$methods(optimize = function() { # no arguments allowed
  "Find posterior mode with the default optimization algorithm"
  stop("FIXME: Implement")
})

StanProgramWithData$methods(sample = function() { # no arguments allowed
  "Sample from posterior distribution with the default MCMC algorithm"
  stop("FIXME: Implement")
})
