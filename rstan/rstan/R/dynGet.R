#  This file is part of RStan
#  Copyright (C) 2015 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

  dynGet <- function (x, 
                      ifnotfound = stop(gettextf("%s not found", sQuote(x)), domain = NA), 
                      minframe = 1L, inherits = FALSE) 
  {
    n <- sys.nframe()
    while (n > minframe) {
      n <- n - 1L
      env <- sys.frame(n)
      if (exists(x, envir = env, inherits = inherits, mode = "numeric")) 
        return(get(x, envir = env, inherits = inherits, mode = "numeric"))
      else if (exists(x, envir = env, inherits = inherits, mode = "logical"))
        return(get(x, envir = env, inherits = inherits, mode = "logical"))
      else if (exists(x, envir = env, inherits = inherits, mode = "list"))
        return(get(x, envir = env, inherits = inherits, mode = "list"))
    }
    return(ifnotfound)
  }
