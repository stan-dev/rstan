CxxFlags <- function(as_character = FALSE) {
  TBB <- system.file("include", .Platform$r_arch, package = "RcppParallel", mustWork = TRUE)
  CXXFLAGS <- paste0("-I, shQuote(TBB), " -D_REENTRANT -DSTAN_THREADS")
  if (isTRUE(as_character)) return(CXXFLAGS)
  cat(CXXFLAGS, " ")
  return(invisible(NULL))
}

LdFlags <- function(as_character = FALSE) {
  TBB <- system.file("lib", .Platform$r_arch, package = "RcppParallel", mustWork = TRUE)
  PKG_LIBS <- paste0("-L", shQuote(TBB), " -Wl,-rpath", TBB, " -ltbb -ltbbmalloc")
  if (isTRUE(as_character)) return(PKG_LIBS)
  cat(PKG_LIBS, " ")
  return(invisible(NULL))
}
