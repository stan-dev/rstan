CXX <- system2(file.path(R.home(component = "bin"), "R"), 
               args = "CMD config CXX", 
               stdout = TRUE, stderr = FALSE)
if (nchar(CXX) > 0) {
  DEV_NULL <- if (.Platform$OS.type == "windows") "nul" else "/dev/null"
  BH <- system.file("include", package = "BH", mustWork = TRUE)
  RcppEigen <- system.file("include", package = "RcppEigen", mustWork = TRUE)
  StanHeaders <- file.path("..", "StanHeaders", "include")
  math <- file.path(StanHeaders, "stan", "math.hpp")
  args <- paste0("-I", BH, " -I", RcppEigen, " -I", StanHeaders, " -o ", DEV_NULL, " ", math)
  check <- system2(CXX, args = args)
  stopifnot(check == 0)
  CXX1X <- system2(file.path(R.home(component = "bin"), "R"), 
                   args = "CMD config CXX1X", 
                   stdout = TRUE, stderr = FALSE)
  check1X <- system2(CXX1X, args = args)
  stopifnot(check1X == 0)
}
