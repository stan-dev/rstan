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

// #include <stan/mcmc/chains.hpp>

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

extern SEXP is_Null_NS(SEXP ns);

#ifdef __cplusplus
}
#endif


/*
 * Tell if it is a NULL native symbol.
 * This function mainly used to tell if a function created by cxxfunction of R
 * package inline points to a NULL address, which would happen when it is
 * deserialized (that is, loaded from what was saved previously by using R's save).
 *
 */
SEXP is_Null_NS(SEXP ns) {
  SEXP ans;
  PROTECT(ans = Rf_allocVector(LGLSXP, 1));
  LOGICAL(ans)[0] = 1;
  PROTECT(ns);
  if (TYPEOF(ns) == EXTPTRSXP) {
    // Rprintf("ptr=%p.\n", EXTPTR_PTR(ns));
    if (EXTPTR_PTR(ns) != NULL) LOGICAL(ans)[0] = 0;
  }
  UNPROTECT(2);
  return ans;
}

