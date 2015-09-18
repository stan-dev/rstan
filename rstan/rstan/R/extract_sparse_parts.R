# This file is part of RStan
# Copyright (C) 2015 Jiqiang Guo, Benjamin Goodrich, and Krzysztof Sakrejda
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

extract_sparse_parts <- function(A, type=c("CSR", "CSC")) {
  type <- match.arg(type)
  if (class(A) == 'matrix') {
    A <- Matrix::Matrix(A, sparse=TRUE)
  }
  if (type == 'CSR') {
    A <- Matrix::t(A)
  }
  if (class(A) == 'dgCMatrix') {
    o <- .Call("extract_sparse_components", A)
  } else {
    stop("Argument 'x' must either be of class 'matrix' or 'dgCMatrix'.")
  }
  return(o)
}




