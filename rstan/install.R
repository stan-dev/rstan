
install_rstan <- function() {
  on.exit(Sys.unsetenv("R_MAKEVARS_USER"))
  on.exit(Sys.unsetenv("R_MAKEVARS_SITE"), add = TRUE)

  opts <- options(repos = c(
    getOption("repos"), 
    rstan = "http://rstan.org/repo/",
    CRAN = "http://cran.rstudio.com"
  ))
  on.exit(options(opts), add = TRUE)
  try(remove.packages("rstan"), silent = TRUE)
  Sys.setenv(R_MAKEVARS_USER = "foobar")
  Sys.setenv(R_MAKEVARS_SITE = "foobar")
  for (pkg in c("inline", "BH", "RcppEigen")) {
    if (!requireNamespace(pkg)) install.packages(pkg)
  }
  if (!requireNamespace("Rcpp")) install.packages("Rcpp", type = "source")
  library(inline) 
  library(Rcpp)
  src <- ' 
    std::vector<std::string> s; 
    s.push_back("hello");
    s.push_back("world");
    return Rcpp::wrap(s);
  '
  hellofun <- cxxfunction(body = src, includes = '', plugin = 'Rcpp', verbose = FALSE)
  test <- try(hellofun())
  if(inherits(test, "try-error")) stop("hello world failed; ask for help on Rcpp list")

  install.packages("rstan", type = 'source')
  library(rstan)
  set_cppo("fast")
  if (any(grepl("^darwin", R.version$os, ignore.case = TRUE))) {
    cat('\nCC=clang', 'CXX=clang++ -arch x86_64 -ftemplate-depth-256', 
        file = "~/.R/Makevars", sep = "\n", append = TRUE)
  }
}
