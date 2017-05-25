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

setMethod("show", signature = "StanType", definition = function(object) {
  if (length(object@values)) cat(str(object@values))
  else cat(class(object))
})

setMethod("summary", signature = "StanType", definition = function(object) {
  leading <- head(dim(object), n = -2L)
  out <- apply(object@values, MARGIN = leading, FUN = fivenum)
  dim(out) <- c(5L, length(out) %/% 5L)
  return(out)
})

# setMethod("show", signature = "StanFactor", definition = function(object) {
#   if( length(dim(object@values)) == 1) {
#     cat("StanType is not filled in yet")
#     return(invisible(NULL))
#   }
#   cat(object@name, "has", tail(dim(object@values), 1), "total draws")
#   return(invisible(NULL))
# })
# 
# setMethod("summary", signature = "StanFactor", definition = function(object) {
#   if( length(dim(object@values)) == 1) {
#     cat("StanType is not filled in yet")
#     return(invisible(NULL))
#   }
#   values <- factor(c(object@values), levels = 1:length(object@labels), 
#                   levels = levels, ordered = object@ordered)
#   tab <- table(values)
#   return(tab)
# })

setMethod("summary", signature = "cov_matrix", definition = function(object) {
  # do not summarize the upper triangle
  dims <- dim(object@values)
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
    out <- rbind(  sapply(1:rows, FUN = monitor_covariance, x = object@values),
                 if(rows > 1) sapply(2:rows, FUN = function(r) {
                   sapply(1:(r-1L), monitor_covariance, r = r, x = object@values)
                 }) )
  }
  else {
    out <- apply(object@values, MARGIN = 1:(len-3L), FUN = function(y) {
      sapply(1:rows, FUN = monitor_covariance, x = y)
    })
    if(rows > 1) out <- rbind(out, apply(object@values, MARGIN = 1:(len-3L),
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
setMethod("Ops", signature = signature("ANY", "StanType"),
          definition = function(e1, e2) {
            return(new(class(e2), values = callGeneric(e1, e2@values)))
})
setMethod("Ops", signature = signature("StanType", "ANY"),
          definition = function(e1, e2) {
            return(new(class(e1), values = callGeneric(e1@values, e2)))
})
setMethod("Ops", signature = signature("StanType", "StanType"),
          definition = function(e1, e2) {
            return(new(class(e1), values = callGeneric(e1@values, e2@values)))
})
setMethod("Math", signature = "StanType", 
          definition = function(x) {
            return(new(class(x), values = callGeneric(x@values)))
})
setMethod("log", signature = "StanType", 
          definition = function(x, ...) {
            return(new(class(x), values = log(x@values, ...)))
})
setMethod("trunc", signature = "StanType", 
          definition = function(x, ...) {
            return(new(class(x), values = trunc(x@values, ...)))
})
setMethod("Math2", signature = "StanType",
          definition = function(x, digits) {
            return(new(class(x), arr = callGeneric(x@values, digits))) 
})
setMethod("Summary", signature = "StanType",
          definition = function(x, ..., na.rm = FALSE) {
            return(new(class(x), values = callGeneric(x@values, ..., na.rm)))
})
setMethod("Complex", signature = "StanType", 
          definition = function(z) {
            return(new(class(z), values = callGeneric(z@values))) 
})

# S3 Methods for a StanType
"[.StanType" <- function(x, i, j, ...) {
  if (missing(j)) return(new(class(x), name = rownames(x@values)[i], 
                             values = x@values[i,], type = x@type))
  else return(values = x@values[i,j,...])
}

"[.real" <- function(x, i = 1L, j, ...) {
  if (is.numeric(i) && i != 1) stop("'i' can only be 1")
  if (is.character(i) && i != x@name)
    stop(paste("'i' can only be", x@name))
  
  if (missing(j)) return(x)
  else return(values = x@values[i,j,...])
}

dim.StanType <- function(x) dim(x@values)

as.array.StanType <- function(x, ...) x@values
as.matrix.StanType <- function(x, ...) {
  dims <- dim(x)
  leading <- head(dims, n = -2L)
  rows <- prod(leading)
  cols <- prod(tail(dims, n =  2L))
  mat <- x@values
  dim(mat) <- c(rows, cols)
  combos <- do.call(expand.grid, args = 
                    lapply(leading, function(z) 1:z))
  rownames(mat) <- paste0(x@name, "[",
                          apply(combos, 1, paste, collapse = ","), "]")
  return(mat)
}
as.matrix.real <- function(x, ...) x@values

