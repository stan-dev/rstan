
library(cOde)
## TODO: switch to cOde code?

## takes the gradient of the expression with respect to vars
symbolic_gradient <- function(expr, vars) lapply(vars, D, expr=expr)
symbolic_jacobian <- function(system_expr, vars) lapply(system_expr, symbolic_gradient, vars=vars)

## generates C++ code from the given jacobian which is in a nested
## list of expressions format.
gen_jac_code <- function(jac, var="J") {
    S <- length(jac)
    V <- length(jac[[1]])

    names(jac) <- 0:(S-1)
    jac <- lapply(jac, function(s) { names(s) <- 0:(V-1); s } )
    jac <- lapply(jac, lapply, function(eq) cOde::replaceOperation("^", "pow", as.character(asOneSidedFormula(eq))[-1]) )
    jac <- unlist(jac)

    names(jac) <- paste0(var, "(", sub("\\.", ",", names(jac)), ")")

    jac_str <- paste0(paste(names(jac), jac, sep="="), ";")
    paste0(paste(jac_str, collapse="\n"),"\n")
}

gen_ode_code <- function(ode_rhs, parameters, export=FALSE, ode_name=deparse(substitute(ode_rhs)), pre_include, post_include) {
    symbols <- unique(unlist(lapply(ode_rhs, all.vars)))
    states <- gsub("^d", "", names(ode_rhs))
    constants <- setdiff(symbols, c(states, parameters))

    cat("Constants :", paste(constants, collapse=", "), "\n")
    cat("Parameters:", paste(parameters, collapse=", "), "\n")
    cat("States    :", paste(states, collapse=", "), "\n")

    S <- length(ode_rhs)
    P <- length(parameters)

    ## compile string for the ode itself
    ode_str <- mapply(function(eq, state) paste0("dydt[", state, "] = ", cOde::replaceOperation("^", "pow", as.character(eq)[-1]), ";\n"), ode_rhs, 1:length(ode_rhs))
    ode_rhs_str <- paste(ode_str, collapse="")

    ode_str_cpp <- mapply(function(eq, state) paste0("dydt(", state, ") = ", cOde::replaceOperation("^", "pow", as.character(eq)[-1]), ";\n"), ode_rhs, 0:(length(ode_rhs)-1))
    ode_rhs_str_cpp <- paste(ode_str_cpp, collapse="")

    if(missing(pre_include) | is.null(pre_include))
        pre_include <- ""
    if(missing(post_include) | is.null(post_include))
        post_include <- ""

    if(file.exists(pre_include))
        ode_rhs_str <- paste0(c(readLines(pre_include), "", ode_rhs_str), collapse="\n")
    if(file.exists(post_include))
        ode_rhs_str <- paste0(c(ode_rhs_str, readLines(post_include), ""), collapse="\n")

    ## compile string for declarations in Stan ...
    def_params <- paste0("real ", parameters, " = theta[", 1:P, "];\n")
    def_states <- paste0("real ", states, " = y[", 1:S, "];\n")
    def_dstates <- paste0("real dydt[", S, "];\n")
    ode_defs <- paste(c(def_dstates, def_states, def_params), collapse="")

    ## ... and in C++
    def_cpp_params <- paste0("const double ", parameters, " = theta_[", 0:(P-1), "];\n")
    def_cpp_states <- paste0("const double ", states, " = y[", 0:(S-1), "];\n")

    ## turn into an expression list suitable for taking derivatives
    ode_rhs_expr <- lapply(ode_rhs, function(ode) parse(text=as.character(ode)[-1]))
    ode_jacobian_states <- symbolic_jacobian(ode_rhs_expr, states)
    ode_jacobian_params <- symbolic_jacobian(ode_rhs_expr, parameters)

    ode_jacobian_states_str <- gen_jac_code(ode_jacobian_states, "Jy")
    ode_jacobian_params_str <- gen_jac_code(ode_jacobian_params, "Jtheta")

    ## check if also c++ pre/post includes are defined
    if(file.exists(sub(".stan", ".hpp", pre_include))) {
        ode_rhs_str_cpp <- paste0(c(readLines(sub(".stan", ".hpp", pre_include)), "", ode_rhs_str_cpp), collapse="\n")
    }
    if(file.exists(sub(".stan", ".hpp", post_include))) {
        ode_rhs_str_cpp <- paste0(c(ode_rhs_str_cpp, "", readLines(sub(".stan", ".hpp", post_include))), collapse="\n")
    }

    stan_def <- paste0("// WARNING: THIS FILE IS AUTO-GENERATED!\n",
                       "real[] ", ode_name, "(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {\n",
                       ode_defs, "\n",
                       ode_rhs_str, "\n",
                       "return(dydt);\n",
                       "}\n"
                       )

    gencode <- list(rhs=ode_rhs_str,
                    rhs_cpp=ode_rhs_str_cpp,
                    name=ode_name,
                    states=states,
                    parameters=parameters,
                    constants=constants,
                    defs=ode_defs,
                    stan_def=stan_def,
                    defs_cpp=defs_cpp <- paste(c(def_cpp_states, def_cpp_params), collapse=""),
                    jac_params=ode_jacobian_params_str,
                    jac_states=ode_jacobian_states_str)

    if(export) {
        cat("Exporting ODE model", ode_name,"...\n")
        ##cat(gencode$defs, file=paste0(ode_name, "_defs.stan"))
        ##cat(gencode$rhs, file=paste0(ode_name, "_rhs.stan"))
        cat(gencode$stan_def, file=paste0(ode_name, ".stan"))

        ## ... and C++ code
        cat(gencode$defs_cpp, file=paste0(ode_name, "_defs.hpp"))
        cat(gencode$rhs_cpp, file=paste0(ode_name, "_rhs.hpp"))
        cat(gencode$jac_params, file=paste0(ode_name, "_jac_params.hpp"))
        cat(gencode$jac_states, file=paste0(ode_name, "_jac_states.hpp"))
    }

    invisible(gencode)
}


## takes Stan ODE model which has appropiate ode definitions
## included. We the parse these, generate the ODE definition and
## generates analytic Jacobians which are written into correct
## includes, parses the stan_file to a stan_file_run which has
## includes applied and finally compiles the Stan model with a on-the
## fly generated Stan parallel ODE code
## extra_magic is for need extra C++ needed for Intel C++
stan_ode_model <- function(stan_ode_file, cpp_ode_dir, compile=TRUE, extra_magic="") {

    ## parse from Stan file the ODE definitions
    ode_def <- grep("\\/\\/ode", readLines(stan_ode_file), value=TRUE)

    if(length(ode_def) == 0)
        stop("Did not find any ODE definition.")

    ode_params <- gsub(" *", "", strsplit(strsplit(grep("ode-parameters", ode_def, value=TRUE), ":")[[1]][2], ",")[[1]])

    if(length(ode_params) == 0)
        stop("No ode-parameters defined!")

    ode_eq <- sapply(strsplit(grep("\\/\\/ode:", ode_def, value=TRUE), "//ode: "), "[", 2)
    ode_eq <- gsub("/dt *=", "=~", ode_eq)
    if(length(ode_eq) == 0)
        stop("No ode equations found!")

    ode_eq <- eval(parse(text=paste0("list(", paste(ode_eq, collapse=", "), ")")))

    extract_filename <- function(str) gsub("(.*: *)([^ ]*)(.*)", "\\2", str)
    pre_include <- extract_filename(grep("\\/\\/ode-pre-include", ode_def, value=TRUE))
    post_include <- extract_filename(grep("\\/\\/ode-post-include", ode_def, value=TRUE))

    if(length(pre_include) == 0)
        pre_include <- NULL
    else if(pre_include == "")
        pre_include <- NULL
    if(length(post_include) == 0)
        post_include <- NULL
    else if(post_include == "")
        post_include <- NULL

    ode_name <- sub(".stan", "_ode", stan_ode_file)
    ode_def_file <- paste0(ode_name, ".stan")
    ode_code <- gen_ode_code(ode_eq, ode_params, export=TRUE, ode_name=ode_name, pre_include=pre_include, post_include=post_include)

    ## paste the parallel integration code into the ode definition file
    ## hmm nested includes seem not to work...
    ##parallel_code <- readLines(file.path(cpp_ode_dir, "integrate_parallel.stan"))
    ##cat(paste0("\n#include \"", file.path(cpp_ode_dir, "integrate_parallel.stan"), "\"\n\n"), file=paste0(ode_name, ".stan"), append=TRUE)
    ##cat(paste0(parallel_code, "\n"), file=ode_def_file, append=TRUE)

    ##ode_name <- deparse(substitute(ode_rhs))
    ##ode_code <- gen_ode_code(ode_rhs, ode_constants, export=TRUE, ode_name=ode_name)

    code <- stanc_builder(file = stan_ode_file, allow_undefined = TRUE)
    ##code$cppcode <- paste0(code$cppcode, extracode)
    ##cat(code$cppcode, file="ode_2cmt_autopar.hpp")

    stan_generated_file <- sub(".stan", "_generated.stan", stan_ode_file)
    cat("Creating", stan_generated_file, "...\n")
    cat(code$model_code, file=stan_generated_file)
    cat("\n", file=stan_generated_file, append=TRUE)

    ode_name <- ode_code$name
    ## set model namespace in line with cmdstan defaults
    ##model_namespace <- sub(".stan", "_model_namespace", stan_generated_file)
    model_namespace <- paste0("model_", sub(".stan", "_namespace", stan_generated_file))
    model_name <- sub(".stan", "", stan_generated_file)

    cpp_extra_include <- paste0(sub(".stan", "", stan_generated_file), "_jacobian_analytic.hpp")

    cpp_defs <- paste0("
namespace model_ns = ", model_namespace, ";
typedef model_ns::", ode_name, "_functor__ ode_functor;
#define ODE_HEADER(header) <", file.path(getwd(), ode_name), "_##header.hpp>
#define MODEL_NAMESPACE ", model_namespace, "
")

    cpp_model_extra <- c(cpp_defs,
                         extra_magic,
                         readLines(file.path(cpp_ode_dir, "coupled_ode_system_jacobian_optional.hpp"))
                         ##readLines(file.path(cpp_ode_dir, "jacobian_ode_analytic.hpp"))
                         )
    cpp_model_extra <- paste(cpp_model_extra, collapse="\n")

    cat("Creating", cpp_extra_include, "...\n")
    cat(cpp_model_extra, file=cpp_extra_include)


#    extracode <- paste0("
#} // escape model namespace
#include \"",  file.path(getwd(), cpp_extra_include), "\"
#namespace ", model_namespace, " {
#")

    if(compile) {
        cat("Compiling ODE Stan model...\n")
        scb <- stanc_builder(stan_generated_file, obfuscate_model_name=FALSE, allow_undefined=TRUE)
        scb$cppcode <- paste0(
            "\n#define STAN_MATH_REV_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_HPP\n",
            scb$cppcode,
            "\n#include \"", file.path(getwd(), cpp_extra_include), "\"\n"
        )
        sm <- stan_model(stanc_ret=scb, model_name=model_name,
                         allow_undefined=TRUE, obfuscate_model_name=FALSE, verbose=FALSE)
        return(sm)
    }

    return(ode_code)
}


