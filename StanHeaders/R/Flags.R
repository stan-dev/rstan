CxxFlags <- function(as_character = FALSE) {
  if (dir.exists(Sys.getenv("TBB_INC"))) {
    TBB_INC <- normalizePath(Sys.getenv("TBB_INC"))
  } else {
    TBB_INC <- system.file("include", package = "RcppParallel", mustWork = TRUE)
  }

  if (file.exists(file.path(TBB_INC, "tbb", "version.h"))) {
    CXXFLAGS <- paste0("-I", shQuote(TBB_INC), " -D_REENTRANT -DSTAN_THREADS -DTBB_INTERFACE_NEW")
  } else {
    CXXFLAGS <- paste0("-I", shQuote(TBB_INC), " -D_REENTRANT -DSTAN_THREADS")
  }

  if (isTRUE(as_character)) return(CXXFLAGS)
  cat(CXXFLAGS, " ")
  return(invisible(NULL))
}

LdFlags <- function(as_character = FALSE) {
  TBB_LIB <- Sys.getenv("TBB_LINK_LIB", Sys.getenv("TBB_LIB"))
  if (dir.exists(TBB_LIB)) {
    TBB_LIB <- normalizePath(TBB_LIB)
  } else {
    TBB_LIB <- system.file("lib", .Platform$r_arch, package = "RcppParallel", mustWork = TRUE)
  }

  PKG_LIBS <- paste0("-L", shQuote(TBB_LIB), " -Wl,-rpath,", shQuote(TBB_LIB), " -ltbb -ltbbmalloc")

  if (isTRUE(as_character)) return(PKG_LIBS)
  cat(PKG_LIBS, " ")
  return(invisible(NULL))
}
