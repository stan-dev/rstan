# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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
                      paste(sort(unique(rosetta$ReturnType)), collapse = ", ")))
  }
  if(is.function(FUN)) FUN <- deparse(substitute(FUN))
  if(!is.character(FUN)) stop("'FUN' must be a character string for a function")
  if(length(FUN) != 1) stop("'FUN' must be of length one")

  if(FUN == "nrow") FUN <- "NROW"
  if(FUN == "ncol") FUN <- "NCOL"

  keep_cols <- colnames(rosetta) != "RFunction"
  if(exists(FUN)) {
    matches <- as.logical(charmatch(rosetta$RFunction, FUN, nomatch = 0L))
    if(any(matches)) return(rosetta[matches, keep_cols, drop = FALSE])
  }
  matches <- grepl(FUN, rosetta$StanFunction, ignore.case = TRUE)
  if(any(matches)) return(rosetta[matches, keep_cols, drop = FALSE])
  else return("no matching Stan functions")
}
