CxxFlags <- function(as_character = FALSE) {
  TBB <- system.file("include", package = "RcppParallel", mustWork = TRUE)
  CXXFLAGS <- paste0("-I", shQuote(TBB), " -D_REENTRANT -DSTAN_THREADS")
  if (isTRUE(as_character)) return(CXXFLAGS)
  cat(CXXFLAGS, " ")
  return(invisible(NULL))
}

LdFlags <- function(as_character = FALSE) {

  # if RcppParallel provides the requisite API for querying
  # the TBB library path, use it
  TBB <- NULL
  if (requireNamespace("RcppParallel", quietly = TRUE)) {
    RcppParallel <- asNamespace("RcppParallel")
    if (is.function(RcppParallel$tbbLibraryPath))
      TBB <- RcppParallel$tbbLibraryPath()
  }

  # Otherwise, find it in the default location (for older RcppParallel releases)
  if (is.null(TBB))
    TBB <- system.file("lib", .Platform$r_arch, package = "RcppParallel", mustWork = TRUE)
  
  PKG_LIBS <- paste0("-L", shQuote(TBB), " -Wl,-rpath,", shQuote(TBB), " -ltbb -ltbbmalloc")
  if (isTRUE(as_character))
    return(PKG_LIBS)

  cat(PKG_LIBS, " ")
  return(invisible(NULL))

}
