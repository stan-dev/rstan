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
