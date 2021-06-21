exposeStanClass <- function(stanc_ret,
                            fields = character(),
                            field_access = c("none", "read_only", "read_write"),
                            where = tempdir(),
                            compile = interactive(), 
                            ...) {
  field_access <- match.arg(field_access)
  if (field_access != "none" && length(fields) == 0L) {
    public <- grep("^public:", stanc_ret$cppcode)
    if (length(public) == 1L) public <- c(grep("^private:", stanc_ret$cppcode), public)
    if (diff(public) > 1L) {
      public[1] <- public[1] + 1L
      public[2] <- public[2] - 1L
      fields <- sub("^.* (.*);$", "\\1", stanc_ret$cppcode[public[1]:public[2]])
    } else fields <- NULL
  }
  if (identical("src", where)) {
    tf <- paste0(stanc_ret$model_cppname, "Module.cc")
    header <- paste0(stanc_ret$model_cppname, ".hpp")
  } else {
    tf <- tempfile(tmpdir = where, fileext = "Module.cc")
    header <- sub("Module\\.cc$", ".hpp", tf)
  }
  
  if (any(grepl("ctor_body", stanc_ret$cppcode))) {
    ctor <- list("rstan::io::rlist_ref_var_context")
  } else ctor <- list(c("SEXP", "SEXP", "SEXP"))
  
  writeLines(stanc_ret$cppcode, con = header)
  meth <- names(stanc_ret$methods)
  Rmeth <- meth
  names(Rmeth) <- sub("_$", "", names(stanc_ret$methods))
  
  FQ <- options()$useFancyQuotes
  options(useFancyQuotes = FALSE)
  on.exit(options(useFancyQuotes = FQ))
  
  Rcpp::exposeClass(class = stanc_ret$model_name,
                    constructors = ctor,
                    fields = fields,
                    methods = meth,
                    file = tf,
                    header = c("// [[Rcpp::depends(rstan)]]",
                               "// [[Rcpp::plugins(cpp14)]]",
                               "#include <RcppEigen.h>",
                               paste0('#include "', header, '"')),
                    CppClass = "stan_model",
                    readOnly = if (field_access == "read_only") fields else character(),
                    rename = Rmeth,
                    Rfile = FALSE)
  if (compile) {
    Rcpp::sourceCpp(file = tf, ...)
    # Rcpp::loadRcppClass(Class = stanc_ret$model_name,
    #                     module = paste0("class_", stanc_ret$model_name),
    #                     where = if (interactive()) .GlobalEnv)
    g <- get(stanc_ret$model_name, envir = .GlobalEnv)@generator
    return(g)
  }
  return(invisible(NULL))
}

