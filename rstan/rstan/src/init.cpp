#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern "C"  SEXP effective_sample_size(SEXP sim, SEXP n_); 
extern "C"  SEXP effective_sample_size2(SEXP sims);
extern "C"  SEXP split_potential_scale_reduction(SEXP sim, SEXP n_); 
extern "C"  SEXP split_potential_scale_reduction2(SEXP sims_);
extern "C"  SEXP seq_permutation(SEXP conf);  
extern "C"  SEXP read_comments(SEXP file, SEXP n);
extern "C"  SEXP stan_prob_autocovariance(SEXP v);
extern "C"  SEXP is_Null_NS(SEXP ns);
extern "C"  SEXP stanc(SEXP model_stancode, SEXP model_name);
extern "C"  SEXP stan_version(); 

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
  CALLDEF(is_Null_NS, 3),
  CALLDEF(effective_sample_size, 2),
  CALLDEF(effective_sample_size2, 1),
  CALLDEF(split_potential_scale_reduction, 2),
  CALLDEF(split_potential_scale_reduction2, 1),
  CALLDEF(seq_permutation, 1),
  CALLDEF(read_comments, 2),
  CALLDEF(stan_prob_autocovariance, 1),
  CALLDEF(stanc, 2),
  CALLDEF(stan_version, 0),
  {NULL, NULL, 0}
};

void attribute_visible R_init_rstan(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
