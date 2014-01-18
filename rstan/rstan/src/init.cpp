/*
 * To register the functions implemented in C++, see 
 * http://cran.r-project.org/doc/manuals/R-exts.html#Registering-native-routines
 *
 * But it seems not to work as it is supposed to be in that
 * even they are still working if not registered, which 
 * is not understood. 
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rversion.h>

#ifdef __cplusplus
extern "C"  {
#endif
SEXP effective_sample_size(SEXP sim, SEXP n_); 
SEXP effective_sample_size2(SEXP sims);
SEXP split_potential_scale_reduction(SEXP sim, SEXP n_); 
SEXP split_potential_scale_reduction2(SEXP sims_);
SEXP seq_permutation(SEXP conf);  
SEXP CPP_read_comments(SEXP file, SEXP n);
SEXP stan_prob_autocovariance(SEXP v);
SEXP is_Null_NS(SEXP ns);
SEXP CPP_stanc(SEXP model_stancode, SEXP model_name);
SEXP CPP_stan_version(); 
#ifdef __cplusplus
}
#endif

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
  CALLDEF(effective_sample_size, 2),
  CALLDEF(effective_sample_size2, 1),
  CALLDEF(split_potential_scale_reduction, 2),
  CALLDEF(split_potential_scale_reduction2, 1),
  CALLDEF(seq_permutation, 1),
  CALLDEF(CPP_read_comments, 2),
  CALLDEF(stan_prob_autocovariance, 1),
  CALLDEF(is_Null_NS, 1),
  CALLDEF(CPP_stanc, 2),
  CALLDEF(CPP_stan_version, 0),
  {NULL, NULL, 0}
};

#ifdef __cplusplus
extern "C"  {
#endif
void attribute_visible R_init_rstan(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  // The call to R_useDynamicSymbols indicates that if the correct C 
  // entry point is not found in the shared library, then an error
  // should be signaled.  Currently, the default
  // behavior in R is to search all other loaded shared libraries for the
  // symbol, which is fairly dangerous behavior.  If you have registered 
  // all routines in your library, then you should set this to FALSE
  // as done in the stats package. [copied from `R Programming for
  // Bioinformatics' // by Robert Gentleman]
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
  R_forceSymbols(dll, TRUE); // copied from package stats, don't know what it does.
#endif
}
#ifdef __cplusplus
}
#endif
