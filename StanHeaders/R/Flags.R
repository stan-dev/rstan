CxxFlags <- function(as_character = FALSE) {
  if (dir.exists(Sys.getenv("TBB_INC"))) {
    TBB_INC <- normalizePath(Sys.getenv("TBB_INC"))
  } else {
    TBB_INC <- system.file("include", package = "RcppParallel", mustWork = TRUE)
  }

  CXXFLAGS <- paste0("-I", shQuote(TBB_INC), " -D_REENTRANT -DSTAN_THREADS")

  if (isTRUE(as_character)) return(CXXFLAGS)
  cat(CXXFLAGS, " ")
  return(invisible(NULL))
}

LdFlags <- function(as_character = FALSE) {
  if (dir.exists(Sys.getenv("TBB_LIB"))) {
    TBB_LIB <- normalizePath(Sys.getenv("TBB_LIB"))
  } else {
    TBB_LIB <- system.file("lib", .Platform$r_arch, package = "RcppParallel", mustWork = TRUE)
  }

  PKG_LIBS <- paste0("-L", shQuote(TBB_LIB), " -Wl,-rpath,", TBB_LIB, " -ltbb -ltbbmalloc")

  if (isTRUE(as_character)) return(PKG_LIBS)
  cat(PKG_LIBS, " ")
  return(invisible(NULL))
}
