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

lookup <- function(FUN, ReturnType = character()) {
  if(length(ReturnType) > 0) {
    if(ReturnType == "ANY") return(rosetta[,-1])
    if(ReturnType %in% rosetta$ReturnType) {
      return(rosetta[rosetta$ReturnType == ReturnType,-1,drop=FALSE])
    }
    else return(paste("no matching Stan functions; ReturnType must be 'ANY' or one of:",
                      paste(sort(unique(rosetta$ReturnType)), collapse = ",")))
  }
  if(is.function(FUN)) FUN <- deparse(substitute(FUN))
  if(!is.character(FUN)) stop("'FUN' must be a character string for a function")
  if(length(FUN) != 1) stop("'FUN' must be of length one")
  
  if(!exists(FUN)) stop(paste("there is no R function by the name of", FUN))
  if(FUN == "nrow") FUN <- "NROW"
  if(FUN == "ncol") FUN <- "NCOL"
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
  
  matches <- as.logical(charmatch(rosetta$RFunction, FUN, nomatch = 0L))
  if(any(matches)) return(rosetta[matches,-1,drop=FALSE])
  matches <- grepl(FUN, rosetta$StanFunction)
  if(any(matches)) return(rosetta[matches,-1,drop=FALSE])
  else return("no matching Stan functions")
}
