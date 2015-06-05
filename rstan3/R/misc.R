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

# a $help method is technically unnecessary but it is more convenient for
# users to ask for help from an instance rather than its ReferenceClass definition
help_from_instance <- function(topic) {
    "Returns documentation string for the $method given by 'topic'"
    CLASS <- .self$getRefClass()
    if (missing(topic)) {
      print(CLASS$help())
    }
    else {
      if (is.name(substitute(topic))) topic <- substitute(topic)
      print(CLASS$help(parse(text = topic)))
    }
    return(invisible(NULL))
}
