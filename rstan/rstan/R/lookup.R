# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015 Jiqiang Guo and Benjamin Goodrich
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

lookup <- function(FUN) {
  if(is.function(FUN)) FUN <- deparse(substitute(FUN))
  if(!is.character(FUN)) stop("'FUN' must be a character string for a function")
  if(length(FUN) != 1) stop("'FUN' must be of length one")
  
  if(FUN == "~") {
    statements <- paste(rosetta$StanFunction[rosetta$Arguments == "~"], "log", sep = "_")
    rosetta <- rosetta[rosetta$StanFunction %in% statements & rosetta$Arguments != "~",]
    rosetta$ReturnType <- NULL    
    colnames(rosetta) <- c("FirstArgument", "StanStatement", "Arguments", "Page")
    rosetta[,1] <- sapply(rosetta$Arguments, FUN = function(x) {
      arg <- scan(text = x, what = character(), sep = ",", quiet = TRUE)[1]
      arg <- sub("(", "", arg, fixed = TRUE)
      return(paste(arg, "~"))
    })
    rosetta$StanStatement <- gsub("_log$", "", rosetta$StanStatement)
    rosetta$Arguments <- gsub("^\\([^,]+, ", "\\(", rosetta$Arguments)
    return(rosetta)
  }
  
  matches <- as.logical(pmatch(rosetta$RFunction, FUN, nomatch = 0L, duplicates.ok = TRUE))
  if(any(matches)) return(rosetta[matches,-1])
  else return("no matching Stan functions")  
}
