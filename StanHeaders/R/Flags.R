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
    # RcppParallel TBB does not have assembly code for Windows ARM64
    # so we need to use compiler builtins
    if (.Platform$OS.type == "windows" && R.version$arch == "aarch64") {
      CXXFLAGS <- paste(CXXFLAGS, "-DTBB_USE_GCC_BUILTINS")
    }
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

  PKG_LIBS <- paste0("-L", shQuote(TBB_LIB))
  # RTools aarch64 does not support rpath, but it is not used on Windows anyway
  if (!(.Platform$OS.type == "windows" && R.version$arch == "aarch64")) {
    PKG_LIBS <- paste0(PKG_LIBS, " -Wl,-rpath,", shQuote(TBB_LIB))
  }
  PKG_LIBS <- paste0(PKG_LIBS, " -ltbb -ltbbmalloc")

  if (isTRUE(as_character)) return(PKG_LIBS)
  cat(PKG_LIBS, " ")
  return(invisible(NULL))
}
