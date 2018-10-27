# This file is part of RStan
# Copyright (C) 2017 Trustees of Columbia University
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

setAs("numeric", "real", function(from) {
  new("real", values = from)
})
setAs("list", "real", function(from) {
  from <- simplify2array(from)
  dims <- dim(from)
  len <- length(dims)
  if (len > 1L) 
    from <- aperm(from, perm = c(dims[len], dims[-len]))
  as(from, "real")
})
setAs("matrix", "real", function(from) {
  new("real", values = as.array(from))
})

setAs("factor", "int", function(from) {
  new("int", values = as.integer(from), 
      labels = levels(from), ordered = is.ordered(from))
})
setAs("numeric", "int", function(from) {
  new("int", values = as.integer(from))
})
setAs("integer", "int", function(from) {
  new("int", values = from)
})
setAs("array", "int", function(from) {
  new("int", values = from)
})

setAs("list", "int", function(from) {
  from <- simplify2array(from)
  dims <- dim(from)
  len <- length(dims)
  if (len > 1L) 
    from <- aperm(from, perm = c(dims[len], dims[-len]))
  as(from, "int")
})

setAs("matrix", "mat", function(from) {
  new("mat", values = from)  
})
setAs("list", "mat", function(from) {
  from <- simplify2array(from)
  dims <- dim(from)
  len <- length(dims)
  # FIXME matrices with 1 row or column
  from <- aperm(from, perm = c(dims[len], dims[-len]))
  as(from, "mat")
})

setAs("numeric", "vec", function(from) {
  dims <- dim(from)
  if (is.null(dims)) dims <- length(from)
  from <- array(from, dim = dims)
  new("vec", values = from)
})
setAs("list", "vec", function(from) {
  from <- simplify2array(from)
  dims <- dim(from)
  len <- length(dims)
  if (len > 1L) 
    from <- aperm(from, perm = c(dims[len], dims[-len]))
  as(from, "vec")
})
