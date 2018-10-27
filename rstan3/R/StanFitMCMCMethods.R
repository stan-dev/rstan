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

StanFitMCMC$methods(show = function() {
  "Show a brief version of the posterior"
  print("brief information about the fit")
  # dims <- dim(.self)
  # cat(dims[1], "unknowns,", dims[2], "chains, and", 
  #     dims[3], "retained samples")
  # return(invisible(NULL))
})

StanFitMCMC$methods(summary = function() {
  "Compute the summary of the posterior"
  print("detailed information about the fit")
  # out <- sapply(.self$sample_draws, FUN = stats4::summary, 
  #               simplify = FALSE)
  # return(t(do.call(cbind, args = out)))
})

StanFitMCMC$methods(as.mcmc.list = function() {
  "Convert to an mcmc.list object for use with the coda package"
  arr <- as.array(.self)
  return(coda::as.mcmc.list(lapply(1:ncol(arr), FUN = function(j) {
    mcmc(t(arr[,j,]))
  })))
})

# S3 methods for StanFitMCMC
"[.StanFitMCMC" <- function(x, i, ...) {
  if (length(i) == 1) {
    out <- x$sample_draws[[i]]
  }
  else {
    out <- x$sample_draws[i]
  }
  return(out)
}

"[[.StanFitMCMC" <- function(x, i, ...) {
  if (is.character(i)) {
    return(sapply(i, simplify = length(i) == 1, FUN = function(y) {
      param <- grep(paste0("^", y, "["), names(x@sample_params))
      return(x@sample_params[[param]][i,])
    }))
  }
  else if (is.numeric(i)) {
    breaks <- cumsum(sapply(x@sample_params, nrow))
    return(sapply(i, simplify = length(i) == 1, FUN = function(y) {
      param <- which(i <= breaks)[1]
      return(x@sample_params[[param]][i - param + 1,])
    }))
  }
}

dim.StanFitMCMC <- function(x) {
  dims <- sapply(x$sample_draws, FUN = dim, simplify = FALSE)
  params <- sum(unlist(sapply(dims, FUN = head, n = -2L)))
  chains <- tail(dims[[1]], 2)
  iterations <- chains[2]
  chains <- chains[1]
  return(c(params, chains, iterations))
}

as.array.StanFitMCMC <- function(x) {
  out <- do.call(rbind, args = 
                   sapply(x$sample_draws, FUN = as.matrix, simplify = FALSE))
  dim(out) <- c(nrow(out), dim(x)[-1])
  return(out)
}

as.matrix.StanFitMCMC <- function(x) {
  out <- as.array(x)
  n <- nrow(out)
  dim(out) <- c(n, length(out) %/% n)
  return(out)
}

names.StanFitMCMC <- function(x) names(x$sample_draws)
