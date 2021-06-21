// This file is part of RStan
// Copyright (C) 2020 Trustees of Columbia University
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

#include <stan/version.hpp>
#include <string>
#include <Rcpp.h>

RcppExport SEXP CPP_stan_version();

SEXP CPP_stan_version() {
  BEGIN_RCPP
  std::string stan_version
  = stan::MAJOR_VERSION + "." +
    stan::MINOR_VERSION + "." +
    stan::PATCH_VERSION;
  SEXP __sexp_result;
  PROTECT(__sexp_result = Rcpp::wrap(stan_version));
  UNPROTECT(1);
  return __sexp_result;
  END_RCPP
}
