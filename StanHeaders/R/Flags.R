CxxFlags <- function(as_character = FALSE) {
  CXXFLAGS <- "-D_REENTRANT -DSTAN_THREADS"
  if (isTRUE(as_character)) return(CXXFLAGS)
  cat(CXXFLAGS, " ")
  return(invisible(NULL))
}

LdFlags <- function(as_character = FALSE) {
  TBB <- system.file("lib", .Platform$r_arch, package = "RcppParallel", mustWork = TRUE)
  if (.Platform$OS.type == "windows") {
    SH <- system.file("libs", .Platform$r_arch, package = "StanHeaders", mustWork = TRUE)
  } else {
    SH <- system.file("lib", .Platform$r_arch, package = "StanHeaders", mustWork = TRUE)
  }
  PKG_LIBS <- paste(paste0("-L", shQuote(TBB)), "-ltbb",
                    paste0("-L", shQuote(SH)), "-lStanHeaders")
  if (isTRUE(as_character)) return(PKG_LIBS)
  cat(PKG_LIBS, " ")
  return(invisible(NULL))
}
