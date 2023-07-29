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

extract_sparse_parts <- function(A) {
  if (!requireNamespace("Matrix")) 
    stop("You have to install the Matrix package to call 'extract_sparse_parts'")
  if (!is(A, 'Matrix')) 
    A <- Matrix::Matrix(A, sparse=TRUE, doDiag=FALSE)
  A <- Matrix::t(A)
  A <- as(A, "dMatrix")
  return(.Call(extract_sparse_components, A))
}
