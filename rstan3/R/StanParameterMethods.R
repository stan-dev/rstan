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

setMethod("initialize", signature = "StanParameter", definition = function(.Object) {
  dims <- dim(.Object@theta)
  for (d in 1:dims) if(dims[d] > 0) dimnames(.Object@theta)[[d]] <- 1:dims[d]
  return(.Object)
})
setMethod("show", signature = "StanParameter", definition = function(object) {
  if( length(dim(object@theta)) == 1) 
    cat("StanParameter is not filled in yet")
  else
    cat("Dimensions of", object@name, "=", head(dim(object@theta), -1), 
        "and has", tail(dim(object@theta), 1), "total draws")
  return(invisible(NULL))
})

setMethod("summary", signature = "StanParameter", definition = function(object) {
  stop("FIXME: Implement")
  monitor(object@theta)
})

setMethod("show", signature = "StanFactor", definition = function(object) {
  if( length(dim(object@theta)) == 1) {
    cat("StanParameter is not filled in yet")
    return(invisible(NULL))
  }
  cat(object@name, "has", tail(dim(object@theta), 1), "total draws")
  return(invisible(NULL))
})

setMethod("summary", signature = "StanFactor", definition = function(object) {
  if( length(dim(object@theta)) == 1) {
    cat("StanParameter is not filled in yet")
    return(invisible(NULL))
  }
  theta <- factor(c(object@theta), levels = 1:length(object@labels), 
                  levels = levels, ordered = object@ordered)
  tab <- table(theta)
  return(tab)
})

setMethod("summary", signature = "StanCovMatrix", definition = function(object) {
  # do not summarize the upper triangle
  dims <- dim(object@theta)
  len <- length(dims)
  rows <- dims[len - 2L]
  cols <- dims[len - 1L]
  iterations <- dims[len]
  chains <- object@chains
  iterations_per_chain <- iterations %/% chains
  monitor_covariance <- function(r, c = r, x) {
    x <- aperm(x, perm = c(iterations, dims[-len]))
    dim(x) <- c(iterations_per_chain, chains, dim(x)[-1L])
    stop("FIXME: Implement monitor")
    monitor(x[,,r,c])
  }
  if (len == 3) {
    out <- rbind(  sapply(1:rows, FUN = monitor_covariance, x = object@theta),
                 if(rows > 1) sapply(2:rows, FUN = function(r) {
                   sapply(1:(r-1L), monitor_covariance, r = r, x = object@theta)
                 }) )
  }
  else {
    out <- apply(object@theta, MARGIN = 1:(len-3L), FUN = function(y) {
      sapply(1:rows, FUN = monitor_covariance, x = y)
    })
    if(rows > 1) out <- rbind(out, apply(object@theta, MARGIN = 1:(len-3L),
                                         FUN = function(y) {
                                sapply(2:rows, FUN = function(r) {
                                  sapply(1:(r - 1L), monitor_covariance(r = r,
                                                                        x = y))
                                })
                              }))
  }
  return(t(out))
})

# (Re)Defining these S4groupGeneric functions facilitates functional programming
setMethod("Ops", signature = signature("ANY", "StanParameter"),
          definition = function(e1, e2) {
            return(new(class(e2), theta = callGeneric(e1, e2@theta)))
})
setMethod("Ops", signature = signature("StanParameter", "ANY"),
          definition = function(e1, e2) {
            return(new(class(e1), theta = callGeneric(e1@theta, e2)))
})
setMethod("Ops", signature = signature("StanParameter", "StanParameter"),
          definition = function(e1, e2) {
            return(new(class(e1), theta = callGeneric(e1@theta, e2@theta)))
})
setMethod("Math", signature = "StanParameter", 
          definition = function(x) {
            return(new(class(x), theta = callGeneric(x@theta)))
})
setMethod("log", signature = "StanParameter", 
          definition = function(x, ...) {
            return(new(class(x), theta = log(x@theta, ...)))
})
setMethod("trunc", signature = "StanParameter", 
          definition = function(x, ...) {
            return(new(class(x), theta = trunc(x@theta, ...)))
})
setMethod("Math2", signature = "StanParameter",
          definition = function(x, digits) {
            return(new(class(x), arr = callGeneric(x@theta, digits))) 
})
setMethod("Summary", signature = "StanParameter",
          definition = function(x, ..., na.rm = FALSE) {
            return(new(class(x), theta = callGeneric(x@theta, ..., na.rm)))
})
setMethod("Complex", signature = "StanParameter", 
          definition = function(z) {
            return(new(class(z), theta = callGeneric(z@theta))) 
})

# S3 Methods for a StanParameter
"[.StanParameter" <- function(x, i, j, ...) {
  if (missing(j)) return(new(class(x), name = rownames(x@theta)[i], 
                             theta = x@theta[i,], type = x@type))
  else return(theta = x@theta[i,j,...])
}

"[.StanReal" <- function(x, i = 1L, j, ...) {
  if (is.numeric(i) && i != 1) stop("'i' can only be 1")
  if (is.character(i) && i != x@name)
    stop(paste("'i' can only be", x@name))
  
  if (missing(j)) return(x)
  else return(theta = x@theta[i,j,...])
}

dim.StanParameter <- function(x) dim(x@theta)

as.array.StanParameter <- function(x, ...) x@theta
as.matrix.StanParameter <- function(x, ...) {
  dims <- dim(x)
  leading <- head(dims, n = -2L)
  rows <- prod(leading)
  cols <- prod(tail(dims, n =  2L))
  mat <- x@theta
  dim(mat) <- c(rows, cols)
  combos <- do.call(expand.grid, args = 
                    lapply(leading, function(z) 1:z))
  rownames(mat) <- paste0(x@name, "[",
                          apply(combos, 1, paste, collapse = ","), "]")
  return(mat)
}
as.matrix.StanReal <- function(x, ...) x@theta

