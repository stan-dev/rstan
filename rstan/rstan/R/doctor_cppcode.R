doctor_cppcode <- function(stanc_ret,
                           ..., # unused but forces exact name matching
                           add_license = character(),
                           use_R_PRNG = FALSE,
                           use_Rcout = FALSE,
                           promote_args_to_auto = FALSE,
                           detemplate = FALSE,
                           double_only = detemplate,
                           propto__ = FALSE,
                           jacobian__ = FALSE,
                           make_data_public = FALSE,
                           default_constructor = FALSE,
                           drop_Eigen = FALSE,
                           drop_log_prob = FALSE,
                           drop_model_header = FALSE,
                           add_toplevel_headers = character(),
                           includes_for_user_defined_functions = character(),
                           return_names = FALSE,
                           return_dims = FALSE,
                           add_methods = list(),
                           add_gradient_method = FALSE,
                           add_hessian_method = FALSE,
                           add_hessian_times_gradient_method = FALSE,
                           add_newton_step = FALSE,
                           methods_for_user_defined_functions = FALSE,
                           friends = character()) {
  dots <- list(...)
  if (length(dots) > 0L) {
    print(str(dots))
    stop("found the above mistakenly passed through ...")
  }
  if (!is.character(add_license))
    stop("'add_license' must be a (possibly empty) character vector")
  if (!is.list(add_methods))
    stop("'add_methods' must be a (possibly empty) list")
  if (!is.character(add_toplevel_headers))
    stop("'add_toplevel_headers' must be a character vector, possibly of length zero")
  if (!is.character(includes_for_user_defined_functions))
    stop("'includes_for_user_defined_functions' must be a character vector, possibly of length zero")
  if (!is.character(friends))
    stop("'friends' must be a (possibly empty) character vector")
  
  if (!is.list(stanc_ret) || is.null(stanc_ret$cppcode) || 
      !is.character(stanc_ret$cppcode))
    stop("'stanc_ret' must be a list (like) that is returned by 'stanc'")
  
  valid <- "[ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_]"
  
  lines <- scan(what = character(), sep = "\n", quiet = TRUE, text = stanc_ret$cppcode)
  
  check_logical_scalar_first <- function(x, na_ok = FALSE) {
    end <- if (na_ok) "must be TRUE, FALSE, or NA" else "must be either TRUE or FALSE"
    msg <- paste(deparse(substitute(x)), end)
    if (!is.logical(x) || length(x) != 1L) stop(msg)
    if (!na_ok && is.na(x)) stop(msg)
    x
  }
  
  four_spaces <- "    "
  protected <- list()
  
  if (check_logical_scalar_first(use_R_PRNG)) {
    protected$PRNG <- paste0(four_spaces, 
                             "boost_random_R base_rng__ = boost_random_R();")
  } else { # this could be dangerous to use with parallel chains
    protected$PRNG <- paste0(four_spaces,
                             "mutable boost::ecuyer1988 base_rng__;")
  }
  
  if (check_logical_scalar_first(use_Rcout)) {
    pstream__ <- "&Rcpp::Rcout"
    lines <- gsub("std::ostream* pstream__ = 0", "std::ostream* pstream__ = &Rcpp::Rcout",
                  lines, fixed = TRUE)
    lines <- gsub("std::ostream* pstream__) const {", 
                  "std::ostream* pstream__ = &Rcpp::Rcout) const {",
                  lines, fixed = TRUE)
  } else {
    pstream__ <- "nullptr"
    lines <- gsub("std::ostream* pstream__ = 0", "std::ostream* pstream__ = nullptr",
                  lines, fixed = TRUE)
    lines <- gsub("std::ostream* pstream__) const {", 
                  "std::ostream* pstream__ = nullptr) const {",
                  lines, fixed = TRUE)
  }
  protected$pstream <- paste0(four_spaces, "std::ostream* pstream__ = ",
                              pstream__, ";")
  
  if (check_logical_scalar_first(promote_args_to_auto)) {
    lines <- gsub("typename boost::math::tools::promote_args<.*$", "auto", lines)
    lines <- gsub("typedef typename boost::math::tools::promote_args<(.*)>.* ",
                  "typedef decltype(\\1) ", lines)
    mark <- grep("typedef decltype", lines)
    lines[mark] <- gsub(", ", " + ", lines[mark], fixed = TRUE)
  }
  
  if (check_logical_scalar_first(detemplate)) {
    if (length(friends) > 0L)
      stop("cannot 'detemplate' the C++ code if 'friends' is specified")
    mark <- grep("template <", lines, fixed = TRUE)
    for (m in rev(mark)) # cannot actually detemplate the RNG for some reason
      if (!grepl("template <typename RNG>", lines[m], fixed = TRUE)) lines <- lines[-m]
      
      if (!isTRUE(double_only)) {
        double_only <- TRUE
        message("'double_only' changed to TRUE due to 'detemplate' being TRUE")
      }
      if (!isTRUE(propto__)) {
        propto__ <- TRUE
        message("'propto__' changed to TRUE due to 'detemplate' being TRUE")
      }
      if (is.na(jacobian__)) {
        jacobian__ <- TRUE
        message("'jacobian__' changed to FALSE due to 'detemplate' being TRUE")
      }
  } else {
    lines <- gsub("typename T__>$", "typename T__ = double>", lines)
    lines <- gsub("typename T_>$",  "typename T_ = double>",  lines)
    lines <- gsub("typename T([[:digit:]]+)__", 
                  "typename T\\1__ = double", lines)
    lines <- gsub("typename T_lp__", "T_lp__ = double", lines)
    lines <- gsub("typename T_lp_accum__", "typename T_lp_accum__ = double", lines)
    lines <- gsub("Class RNG", 
                  paste0("Class RNG = ", ifelse(use_R_PRNG, "boost_random_R",
                                                "boost::ecuyer1988")), lines)
  }
  
  if (check_logical_scalar_first(double_only)) {
    if (!detemplate)
      stop("if 'double_only' is TRUE, then 'detemplate' must also be TRUE")
    lines <- gsub("T_lp_accum__", "accumulator<double>", lines)
    lines <- gsub("T_lp__", "double", lines)
    lines <- gsub("T[[:digit:]]*_+", "double", lines)
  }
  
  if (check_logical_scalar_first(propto__) || !propto__) {
    if (!detemplate) {
      lines <- gsub("bool propto__", 
                    paste("bool propto__ =", ifelse(propto__, "true", "false")),
                    lines, fixed = TRUE)
      lines <- gsub("bool propto,", 
                    paste("bool propto =", ifelse(propto__, "true,", "false,")),
                    lines, fixed = TRUE)
      mark <- grep("const static bool propto", lines)
      lines[mark] <- paste0(four_spaces, "const static bool propto__ = true;")
    }
    else {
      protected$propto__ <- paste0(four_spaces, "const bool propto__ = ",
                                   ifelse(propto__, "true", "false"), ";")
      protected$propto <- paste0(four_spaces, "const bool propto = ",
                                 ifelse(propto__, "true", "false"), ";")
    }
  }
  
  if (check_logical_scalar_first(jacobian__) || !jacobian__) {
    if (!detemplate) {
      lines <- gsub("bool jacobian__",
                    paste("bool jacobian__ =", ifelse(jacobian__, "true", "false")),
                    lines, fixed = TRUE)
      lines <- gsub("bool jacobian,",
                    paste("bool jacobian =", ifelse(jacobian__, "true,", "false,")),
                    lines, fixed = TRUE)
    }
    else {
      protected$jacobian__ <- paste0(four_spaces, "const bool jacobian__ = ",
                                     ifelse(jacobian__, "true", "false"), ";")
      protected$jacobian <- paste0(four_spaces, "const bool jacobian = ",
                                   ifelse(jacobian__, "true", "false"), ";")
    }
  }
  
  if (check_logical_scalar_first(default_constructor)) {
    if (!isTRUE(make_data_public)) {
      make_data_public <- TRUE
      message("'make_data_public' changed to TRUE due to 'default_constructor' being TRUE")
    }
  }
  
  if (check_logical_scalar_first(make_data_public)) {
    lines <- sub("^private:", "public:", lines)
  }
  
  find_closing <- function(lines, mark) {
    if (length(mark) != 1 || is.na(mark)) stop("mark is wrong")
    opening <- lines[mark]
    chars <- strsplit(opening, split = NULL, fixed = TRUE)[[1]]
    count <- 0L
    pos <- 1L
    while (pos <= length(chars) && chars[pos] == " ") {
      count <- count + 1L
      pos <- pos + 1L
    }
    closing <- grep(paste0("^[[:space:]]{", count, "}", "\\}"), lines[-c(1:mark)])[1]
    if (length(closing) == 0) stop("no closing found")
    closing <- closing + mark
    return(closing)
  }
  
  if (check_logical_scalar_first(drop_Eigen)) {
    mark <- grep("log_prob(Eigen::Matrix<", lines, fixed = TRUE)
    if (!detemplate) mark <- mark - 1L
    closing <- find_closing(lines, mark)
    lines <- lines[-c(mark:closing)]
    mark <- grep("Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,", 
                 lines, fixed = TRUE)[1] - 1L
    closing <- find_closing(lines, mark)
    lines <- lines[-c(mark:closing)]
    mark <- grep("Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,", 
                 lines, fixed = TRUE)[1] - 1L
    closing <- find_closing(lines, mark)
    lines <- lines[-c(mark:closing)]
  }
  
  if (check_logical_scalar_first(drop_log_prob)) {
    mark <- grep("log_prob(", lines, fixed = TRUE)
    mark <- mark[!grepl("}", lines[mark])]
    if (!detemplate) mark <- mark - 1L
    closing <- find_closing(lines, mark)
    lines <- lines[-c(mark:closing)]
  }
  
  necessary_headers <- c("#include <rstan/boost_random_R.hpp>",
                         "#include <rstan/io/rlist_ref_var_context.hpp>")
  if (check_logical_scalar_first(drop_model_header)) {
    lines <- grep("#include <stan/model/model_header.hpp>", lines,
                  fixed = TRUE, invert = TRUE, value = TRUE)
    necessary_headers <- c(necessary_headers,
                           ifelse(double_only, 
                                  "#include <stan/math/prim/mat.hpp>",
                                  "#include <stan/math/math.hpp>"),
                           "#include <stan/model/prob_grad.hpp>",
                           "#include <stan/io/array_var_context.hpp>",
                           "#include <stan/io/cmd_line.hpp>",
                           "#include <stan/io/dump.hpp>",
                           "#include <stan/io/reader.hpp>",
                           "#include <stan/io/writer.hpp>",
                           "#include <stan/lang/rethrow_located.hpp>",
                           "#include <stan/model/prob_grad.hpp>",
                           "#include <stan/model/indexing.hpp>",
                           "#include <boost/exception/all.hpp>")
    if (!use_R_PRNG)
      necessary_headers <- c(necessary_headers,
                             "#include <boost/random/additive_combine.hpp>",
                             "#include <boost/random/linear_congruential.hpp>")
    
  }
  
  if (check_logical_scalar_first(return_names)) {
    get_param_names <- c(
      "() const {",
      "  std::vector<std::string> names__;",
      "  get_param_names(names__);",
      "  return names__;",
      "}")
    constrained_param_names <- c(
      "(bool include_tparams__ = true, bool include_gqs__ = true) const {",
      "  std::vector<std::string> param_names__;",
      "  constrained_param_names(param_names__, include_tparams__, include_gqs__);",
      "  return param_names__;",
      "}")
    unconstrained_param_names <- sub("constrained_param_names",
                                     "unconstrained_param_names", 
                                     constrained_param_names)
    add_methods$get_param_names_ <- get_param_names
    add_methods$constrained_param_names_ <- constrained_param_names
    add_methods$unconstrained_param_names_ <- unconstrained_param_names
  }
  
  if (check_logical_scalar_first(return_dims)) {
    get_dims <- c(
      "() const {",
      "  std::vector<std::vector<size_t> > dimss__;",
      "  get_dims(dimss__);",
      "  return dimss__;",
      "}")
    add_methods$get_dims_ <- get_dims
  }
  
  if (!isTRUE(use_R_PRNG)) {
    set_Boost_seed <- c(
      "(const unsigned int seed) const {",
      "  base_rng__.seed(seed);",
      "  return;",
      "}")
    add_methods$set_Boost_seed <- set_Boost_seed
  }
  
  if (!drop_log_prob) { # this is easily exportable to R
    log_prob <- if (detemplate) c(
      "(std::vector<double>& params_r__) const {",
      "  std::vector<int> params_i__;",
      "  return log_prob(params_r__, params_i__, pstream__);",
      "}") else c(
        "(std::vector<double>& params_r__) const {",
        "  std::vector<int> params_i__;",
        "  return log_prob<>(params_r__, params_i__, pstream__);",
        "}")
          
    add_methods$log_prob_ <- log_prob
  }
  
  if (check_logical_scalar_first(add_gradient_method)) {
    if (drop_log_prob)
      stop("if 'add_gradient_method' is TRUE, 'drop_log_prob' must be FALSE")
    if (double_only)
      stop("if 'add_gradient_method' is TRUE, 'double_only' must be FALSE")
    gradient <- c(
      "(Eigen::VectorXd params_r__) {",
      "  std::vector<int> params_i__;",
      "  double fx;",
      "  std::vector<double> grad_fx;",
      "  auto stan::math::gradient([&params_i__](auto theta) {",
      "    return log_prob<propto__, jacobian__>(to_array_1d(theta), params_i__, pstream__);",
      "  }, params_r__, fx, grad_fx);",
      "  return grad_fx;",
      "}")
    add_methods$gradient <- gradient
  }

  if (check_logical_scalar_first(add_hessian_method)) {
    if (drop_log_prob) 
      stop("if 'add_gradient_method' is TRUE, 'drop_log_prob' must be FALSE")
    if (detemplate)
      stop("if 'add_gradient_method' is TRUE, 'detemplate' must be FALSE")
    hessian <- c(
      "(std::vector<double> params_r__) {",
      "  std::vector<int> params_i__;",
      "  std::vector<double> grad_fx;",
      "  std::vector<double> hess_fx;",
      paste0("  auto stan::model::grad_hess_log_prob<propto__, jacobian__>(", 
             "this, params_r__, params_i__, grad_fx, hess_fx, pstream__);"),
      "  return hess_fx;",
      "}")
    add_methods$hessian <- hessian
    necessary_headers <- c(necessary_headers,
                           "#include <stan/model/grad_hess_log_prob.hpp>")
  }
  
  if (check_logical_scalar_first(add_hessian_times_gradient_method)) {
    if (drop_log_prob) 
      stop("if 'add_hessian_times_gradient_method' is TRUE, 'drop_log_prob' must be FALSE")
    if (detemplate)
      stop("if 'add_hessian_times_gradient_method' is TRUE, 'detemplate' must be FALSE")
    hessian_times_gradient <- c(
      "(std::vector<double> params_r__) {",
      "  std::vector<int> params_i__;",
      "  std::vector<double> grad_fx;",
      "  std::vector<double> hess_fx;",
      paste0("  auto stan::model::grad_hess_log_prob<propto__, jacobian__>(", 
             "this, params_r__, params_i__, grad_fx, hess_fx, pstream__);"),
             "  return (hess_fx.selfadjointView<Eigen::Lower>() * grad_fx).eval();",
      "}")
    add_methods$hessian_times_gradient <- hessian_times_gradient
    if (!add_hessian_method)
      necessary_headers <- c(necessary_headers,
                             "#include <stan/model/grad_hess_log_prob.hpp>")
  }
  
  if (check_logical_scalar_first(add_newton_step)) {
    if (isTRUE(drop_log_prob)) 
      stop("if 'add_newton_step' is TRUE, 'drop_log_prob' must be FALSE")
    newton_step <- c(
      "(std::vector<double> params_r__) {",
      "  std::vector<int> params_i__;",
      paste0("  auto stan::optimization::newton_step",
             "(this, params_r, params_i, pstream__);"),
      "  return params_r__;",
      "}")
    add_methods$newton_step <- newton_step
    necessary_headers <- c(necessary_headers, 
                           "#include <stan/services/optimization/newton_step.hpp>")
  }
  
  make_sig_bod <- function(sig, friend) {
    if (detemplate) {
      start <- if (friend) "(this, " else "("
    } else start <- if (friend) "<>(this, " else "<>("
    if (grepl("lp__", sig)) {
      sig <- sub("(^.*lp__).*$", "\\1) {", sig)
      inputs <- strsplit(sig, split = ",")[[1]]
      inputs <- sub(") {", "", inputs, fixed = TRUE)
      inputs <- sub(paste0("^.* (", valid, "+$)"), "\\1", inputs)
      if (friend) inputs <- inputs[-1L]
      inputs <- c(inputs, "lp_accum__", "pstream__") # dubious
      inputs <- paste(inputs, collapse = ", ")
      bod <- c("accumulator<double> lp_accum__;",
               paste0("return ", namespace_name, "::", method_name,
                      start, inputs, "); }"))
    } else if (grepl("base_rng__", sig)) {
      sig <- sub("(^.*), RNG.*$", "\\1) {", sig)
      inputs <- strsplit(sig, split = ",")[[1]]
      inputs <- sub(") {", "", inputs, fixed = TRUE)
      inputs <- sub(paste0("^.* (", valid, "+$)"), "\\1", inputs)
      if (friend) inputs <- inputs[-1L]
      inputs <- c(inputs, "base_rng__", "pstream__")
      inputs <- paste(inputs, collapse = ", ")
      bod <- paste0("return ", namespace_name, "::", method_name,
                    start, inputs, "); }")
    } else if (grepl(",", sig)) {
      sig <- sub("(^.*), std::ostream.*$", "\\1) {", sig)
      inputs <- strsplit(sig, split = ",")[[1]]
      inputs <- sub(") {", "", inputs, fixed = TRUE)
      inputs <- sub(paste0("^.* (", valid, "+$)"), "\\1", inputs)
      if (friend) inputs <- inputs[-1L]
      inputs <- c(inputs, "pstream__")
      inputs <- paste(inputs, collapse = ", ")
      bod <- paste0("return ", namespace_name, "::", method_name,
                    start, inputs, "); }")
    } else {
      if (friend) stop("no template arguments in friend function")
      sig <- "() {"
      bod <- paste0("return ", namespace_name, "::", method_name, "(); }")
    }
    return(list(sig = sig, bod = bod))    
  }
  
  lines <- gsub("Eigen::Matrix<double, Eigen::Dynamic, 1>", 
                "Eigen::VectorXd", lines, fixed = TRUE)
  lines <- gsub("Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>", 
                "Eigen::MatrixXd", lines, fixed = TRUE)

  mark <- grep("^namespace .* \\{$", lines)
  namespace_name <- sub("namespace ", "", lines[mark])
  namespace_name <- sub(" {", "", namespace_name, fixed = TRUE)
  if (check_logical_scalar_first(methods_for_user_defined_functions)) {
    mark <- grep(paste0("^", valid, "+\\("), lines)
    if (length(mark) == 0L)
      stop("'methods_for_user_defined_functions' is TRUE but no functions found")
    for (i in rev(mark)) {
      end <- i
      while (!grepl("{", lines[end], fixed = TRUE)) end <- end + 1L
      sig <- paste(lines[i:end], collapse = " ")
      method_name <- sub(paste0("^(", valid, "+)\\(.*$"), "\\1", sig)
      sig <- sub(method_name, "", sig, fixed = TRUE)
      sig_bod <- make_sig_bod(sig, friend = FALSE)
      if (!detemplate) 
        sig_body$sig <- gsub("T[[:digit:]]+__", "double", sig_bod$sig)
      add_methods[[method_name]] <- c(sig_bod$sig, sig_bod$bod)
    }
  }
  
  if (length(friends) > 0L) {
    mark <- grep(paste0("^", valid, "+\\(.*\\);$"), lines)
    if (length(mark) == 0L)
      stop("'friends' is specified but no undefined functions found")
    friends_list <- list()
    for (i in rev(mark)) {
      sig <- lines[i]
      method_name <- sub(paste0("^(", valid, "+\\(.*$"), "\\1", sig)
      if (!(method_name %*% names(friends))) next
      sig <- sub(method_name, "", sig, fixed = TRUE)
      sig_body <- make_sig_bod(sig, friend = TRUE)
      add_methods[[method_name]] <- c(sig_bod$sig, sig_bod$bod)
      sig <- sub("(", "(prob_grad&, ", sig_bod$sig, fixed = TRUE)
      sig <- sub("{", ";", sig, fixed = TRUE)
      friends_list[[method_name]] <- sig
    }
  }
  
  for (i in seq_along(add_methods)) {
    mark <- grep("^}; // model$", lines) - 1L
    method <- add_methods[[i]]
    lines <- append(lines, after = mark,
                    values = c(paste0(four_spaces, "auto ", 
                                      names(add_methods)[i], method[1]), 
                               paste0(four_spaces, four_spaces, method[-1])))
  }
  
  if (length(friends)) {
    mark <- grep("^}; // model$", lines) - 1L
    lines <- append(lines, values = "friend:", after = mark)
    for (i in seq_along(friends_list)) {
      mark <- grep("^}; // model$", lines) - 1L
      method <- friends_list[[i]]
      lines <- append(lines, after = mark,
                      values = c(paste0(four_spaces, "auto ", 
                                        names(friends_list)[i], method[1]), 
                                 paste0(four_spaces, four_spaces, method[-1])))
    }
  }

  # do not create base_rng__ in ctor_body
  mark <- grep("boost::ecuyer1988 base_rng__ =", lines, fixed = TRUE)
  lines <- lines[-c(mark:(mark + 2L))]
  
  # deal with constructor
  mark <- grep(paste0("^class[[:space:]]+", 
                      "[ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz]", 
                      valid, "*[[:space:]]*: public prob_grad \\{"), lines)
  class_name <- sub("^class (.*) : public prob_grad \\{", "\\1", lines[mark])
  lines <- append(lines, after = mark,
                  values = c("protected:", unlist(protected)))
  lines <- append(lines, after = mark - 1L,
                  values = includes_for_user_defined_functions)
  mark <- grep(paste0(class_name, "(stan::io::var_context& context__,"), 
               lines, fixed = TRUE)[1]
  closing <- find_closing(lines, mark)
  lines <- lines[-c(mark:closing)]
  mark <- grep(paste0(class_name, "(stan::io::var_context& context__,"), 
               lines, fixed = TRUE)[1]
  closing <- find_closing(lines, mark)
  lines <- lines[-c(mark:closing)]
  if (default_constructor) {
    mark <- grep("void ctor_body(stan::io::var_context& context__,", 
                 lines, fixed = TRUE)
    closing <- find_closing(lines, mark)
    lines <- lines[-c(mark:closing)]
    ctor <- paste0(four_spaces, class_name, "() : prob_grad(0) { }")
  } else {
    ctor <- c(paste0(four_spaces, class_name, 
                     "(rstan::io::rlist_ref_var_context context__) ",
                     ": prob_grad(0) {"), 
              paste0(four_spaces, four_spaces,
                     "ctor_body(context__, 0, ", pstream__, ");"),
              paste0(four_spaces, "}"))
  }
  mark <- grep("^public:", lines)[1L + make_data_public]
  lines <- append(lines, values = ctor, after = mark)
  lines <- append(lines, after = 1L, c(necessary_headers, 
                                       add_toplevel_headers))

  if (length(add_license) == 1L) {
    license <- file.path(R.home("share"), "licenses", add_license)
    if (!file.exists(license))
      stop("'add_license' not among license abbreviations known to R")
    lines < c("/*", readLines(license), "*/", "", lines)
  }
  
  stanc_ret$cppcode <- lines
  stanc_ret$methods <- add_methods
  stanc_ret$friends <- if (length(friends)) friends_list
  return(stanc_ret)
}
