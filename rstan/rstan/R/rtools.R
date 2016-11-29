# This file is part of RStan
# Copyright (C) 2015, 2016 Trustees of Columbia University
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

find_rtools <- function(...) {
  if (.Platform$OS.type == "unix") return(TRUE)
  CXX <- Sys.which("g++")
  if (nchar(CXX) > 0) return(TRUE)
  CXX <- system2(file.path(Sys.getenv("R_HOME"), "bin",
                           Sys.getenv("R_ARCH_BIN"), "R"),
                 args = "CMD config CXX", stdout = TRUE)
  if (!is.null(attributes(CXX)$status)) {
    # do nothing
  }
  else if (nchar(CXX) > 0) return(TRUE)
  message("No C++ compiler found, so the following will probably not work.")
  message("See https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows#toolchain")
  return(FALSE)
}
