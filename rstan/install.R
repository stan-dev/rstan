
install_rstan <- function() {
  on.exit(Sys.unsetenv("R_MAKEVARS_USER"))
  on.exit(Sys.unsetenv("R_MAKEVARS_SITE"), add = TRUE)

  try(remove.packages("rstan"), silent = TRUE)
  Sys.setenv(R_MAKEVARS_USER = "foobar")
  Sys.setenv(R_MAKEVARS_SITE = "foobar")
  install.packages(c("inline", "BH", "RcppEigen"))
  install.packages("Rcpp", type = "source")
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

  options(repos = c(getOption("repos"), 
          rstan = "http://rstan.org/repo/"))
  install.packages("rstan", type = 'source')
  library(rstan)
  set_cppo("fast")
  if (any(grepl("^darwin", R.version$os, ignore.case = TRUE))) {
    cat('\nCC=clang', 'CXX=clang++ -arch x86_64 -ftemplate-depth-256', 
        file = "~/.R/Makevars", sep = "\n", append = TRUE)
  }
  return(invisible(NULL))
}
