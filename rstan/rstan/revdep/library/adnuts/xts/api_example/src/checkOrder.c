/*
Example of using the C API from xts in new package code

This is derived from examining the source of 
packages Matrix and lme4.
*/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <R_ext/Rdynload.h>  /* required by R */
/*
  The following header file is from the include directory that is
  included with xts
*/
#include "xtsAPI.h"  /* function declaration and macros */


SEXP check_order (SEXP x, SEXP incr, SEXP strict)
{
  SEXP ret;
  PROTECT(ret = allocVector(LGLSXP, 1));
  /*
     do_is_ordered is imported from the xts package.
     All that is needed here is to call it.
  */
  ret = xtsIsOrdered(x, incr, strict);
  UNPROTECT(1);
  return ret;
}
