CXX <- system2(file.path(R.home(component = "bin"), "R"), 
               args = "CMD config CXX", 
               stdout = TRUE, stderr = FALSE)
if (nchar(CXX) > 0) {
  DEV_NULL <- if (.Platform$OS.type == "windows") "nul" else "/dev/null"
  BH <- system.file("include", package = "BH", mustWork = TRUE)
  RcppEigen <- system.file("include", package = "RcppEigen", mustWork = TRUE)
  StanHeaders <- system.file("include", package = "StanHeaders", mustWork = TRUE)
  math <- file.path(StanHeaders, "stan", "math.hpp")
  m64 <- grepl("m64", CXX)
  if (m64) CXX <- sub(" -m64", "", CXX)
  args <- paste0(if(m64) "-m64 ", "-I", BH, " -I", RcppEigen, " -I", StanHeaders, 
                 " -o ", DEV_NULL, " ", math)
  check <- system2(CXX, args = args)
  if (check != 0) system2(CXX, args = args, stdout = TRUE, stderr = TRUE)
  CXX1X <- system2(file.path(R.home(component = "bin"), "R"), 
                   args = "CMD config CXX1X", 
                   stdout = TRUE, stderr = FALSE)
  if (nchar(CXX1X) > 0) {
    m64 <- grepl("m64", CXX1X)
    if (m64) CXX1X <- sub(" -m64", "", CXX1X)
    CXX1XSTD <- system2(file.path(R.home(component = "bin"), "R"), 
                        args = "CMD config CXX1XSTD", 
                        stdout = TRUE, stderr = FALSE)
    args <- paste0(if(m64) "-m64 ", CXX1XSTD, " -I", BH, " -I", RcppEigen, " -I", StanHeaders,
                   " -o ", DEV_NULL, " ", math)
    check1X <- system2(CXX1X, args = args)
    if (check1X == 0) system2(CXX1X, args = args, stdout = TRUE, stderr = TRUE)
    if (.Platform$OS.type == "windows")
      warning("Compilation is expected to fail under 'g++ -std=gnu++0x'")
  }
}
