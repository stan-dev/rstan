// This file is part of RStan
// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
//
// RStan is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// RStan is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

/*
 * To register the functions implemented in C++, see
 * http://cran.r-project.org/doc/manuals/R-exts.html#Registering-native-routines
 *
 * But it seems not to work as it is supposed to be in that
 * they are still working if not registered.
 */
#include <RcppEigen.h>
#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rversion.h>


using namespace Rcpp;

RcppExport SEXP _rcpp_module_boot_class_model_base();

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
SEXP CPP_stanc280(SEXP model_stancode, SEXP model_name, 
                  SEXP allow_undefined, SEXP include_paths);
SEXP stanfuncs(SEXP model_stancode, SEXP model_name, SEXP allow_undefined);
SEXP CPP_stan_version();
SEXP extract_sparse_components(SEXP A);
SEXP get_rng_(SEXP seed);
SEXP get_stream_();

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
  CALLDEF(CPP_stanc280, 4),
  CALLDEF(stanfuncs, 3),
  CALLDEF(CPP_stan_version, 0),
  CALLDEF(extract_sparse_components, 1),
  CALLDEF(get_rng_, 1),
  CALLDEF(get_stream_, 0),
  {"_rcpp_module_boot_class_model_base", (DL_FUNC) &_rcpp_module_boot_class_model_base, 0},
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
}
#ifdef __cplusplus
}
#endif
