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

rstan.package.skeleton <- function(name = "anRpackage", list = character(),
                                   environment = .GlobalEnv,
                                   path = ".", force = FALSE,
                                   code_files = character(),
                                   stan_files = character()) {
  stop(
    "'rstan.package.skeleton' is now 'rstan_package_skeleton' ",
    "and can be found in the 'rstantools' package.",
    call. = FALSE
  )
}
