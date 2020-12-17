stanFunction <- function(function_name, ..., env = parent.frame(), rebuild = FALSE,
                         cacheDir = getOption("rcpp.cache.dir", tempdir()), 
                         showOutput = verbose, verbose = getOption("verbose")) {
  make_type <- function(x, recursive = FALSE) {
    is_array <- is.list(x)
    if (is_array) {
      base_type <- make_type(x[[1L]], recursive = TRUE)
      if (recursive) return(base_type)
      type <- sub("const ", "", base_type)
      j <- 1L
      while(j <= length(x) && is.list(x[[j]])) {
        type <- paste0("std::vector<", type, " >")
        j <- j + 1L
      }
      type <- paste0("const std::vector<", type, " >&")
      return(type)
    }
    Eigen <- FALSE
    if (is.matrix(x)) {
      Eigen <- TRUE
      if (nrow(x) == 1L) type <- "stan::math::row_vector_d"
      else type <- "stan::math::matrix_d"
    } else if (length(x) > 1L) {
      if (is.integer(x)) {
        type <- "std::vector<int>"
      } else {
        Eigen <- TRUE
        type <- "stan::math::vector_d"
      }
    } else if (is.integer(x)) {
      type <- "int"
    } else if (is.numeric(x)) {
      type <- "double"
    } else stop(paste("all arguments to", function_name, "must be matrices,",
                      "vectors, integers, doubles or lists thereof"))
    if (Eigen) type <- paste0("const ", type, "&")
    else type <- paste0("const ", type)
    return(type)
  }
  DOTS <- list(...)
  types <- sapply(DOTS, FUN = make_type)
  double_lists <- types == "const std::vector<double >&"
  if (any(double_lists)) types[double_lists] <- "const List&"
  int_lists <- types == "const std::vector<int >&"
  if (any(int_lists)) types[int_lists] <- "const List&"
  code <- paste0("auto ", function_name, "(",
                 paste(types, names(types), collapse = ", "), 
                 ") { return stan::math::", function_name, "(", 
                 paste(ifelse(double_lists,
                              paste0("std::vector<double>(", names(types), ".begin(), ",
                                                             names(types), ".end())"),
                              ifelse(int_lists,
                                     paste0("std::vector<int>(", names(types), ".begin(), ",
                                                                 names(types), ".end())"),
                                     names(types))), collapse = ", "), "); }")
  incl <- dir(system.file("include", "stan", "math", "prim", 
                          package = "StanHeaders", mustWork = TRUE),
              pattern = "hpp$")
  incl <- setdiff(incl, "core.hpp")
  incl <- paste0("#include <stan/math/prim/", incl, ">")
  if (grepl("_rng$", function_name)) {
    create_rng <- system.file("include", "src", "stan", "services", "util", "create_rng.hpp",
                              package = "StanHeaders", mustWork = TRUE)
    op <- options("useFancyQuotes")
    options(useFancyQuotes = FALSE)
    on.exit(options(useFancyQuotes = op))
    incl <- c(incl, paste0('#include ', dQuote(create_rng)))
    code <- sub(") {", ", const int random_seed = 0) {", code, fixed = TRUE)
    code <- sub(" return ", 
                "boost::ecuyer1988 base_rng__ = stan::services::util::create_rng(random_seed, 0); return ",
                code)
      code <- sub("); }", ", base_rng__); }", code, fixed = TRUE)
  }
  old_USE_CXX14 <- Sys.getenv("USE_CXX14")
  on.exit(Sys.setenv(USE_CXX14 = old_USE_CXX14))
  Sys.setenv(USE_CXX14 = "1")
  Rcpp::cppFunction(code, depends = c("StanHeaders", "RcppEigen", "BH"),
                    includes = incl, env = env, rebuild = rebuild, 
                    cacheDir = cacheDir,
                    showOutput = showOutput, verbose = verbose)
  if (grepl("_rng$", function_name)) {
    fun <- get(function_name, envir = env, mode = "function")
    formals(fun)$random_seed <- quote(sample.int(.Machine$integer.max, size = 1L))
    assign(function_name, value = fun, envir = env)
  }
  return(do.call(function_name, args = DOTS, envir = env))
}

