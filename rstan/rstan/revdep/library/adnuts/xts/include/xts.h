/*
Header file for using internal C-level facilities
provided by xts.

This is not 100% designed for end users, so
any user comments and bug reports are very
welcomed.

Copyright Jeffrey A. Ryan 2008

This source is distributed with the same license
as the full xts software, GPL (>= 2).
*/
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifndef _XTS
#define _XTS

#ifdef __cplusplus
extern "C" {
#endif

/*
INTERNAL SYMBOLS
*/
SEXP xts_IndexSymbol;
SEXP xts_ClassSymbol;
SEXP xts_IndexFormatSymbol;
SEXP xts_IndexClassSymbol;
SEXP xts_IndexTZSymbol;
SEXP xts_IndexTclassSymbol;
SEXP xts_IndexTzoneSymbol;

/*
DATA TOOLS
*/
#define  xts_ATTRIB(x)                  coerceVector(do_xtsAttributes(x),LISTSXP)
#define  xts_COREATTRIB(x)              coerceVector(do_xtsCoreAttributes(x),LISTSXP)

// attr(x, 'index') or .index(x)
#define  GET_xtsIndex(x)                getAttrib(x, xts_IndexSymbol)
#define  SET_xtsIndex(x,value)          setAttrib(x, xts_IndexSymbol, value)

// attr(x, '.indexCLASS') or indexClass(x)
#define  GET_xtsIndexClass(x)           getAttrib(x, xts_IndexClassSymbol)
#define  SET_xtsIndexClass(x,value)     setAttrib(x, xts_IndexvalueSymbol, value)

// attr(x, '.indexFORMAT') or indexFormat(x)
#define  GET_xtsIndexFormat(x)          getAttrib(x, xts_IndexFormatSymbol)
#define  SET_xtsIndexFormat(x,value)    setAttrib(x, xts_IndexFormatSymbol, value)

// attr(x, '.indexTZ') or indexTZ(x)
#define  GET_xtsIndexTZ(x)              getAttrib(x, xts_IndexTZSymbol)
#define  SET_xtsIndexTZ(x,value)        setAttrib(x, xts_IndexTZSymbol, value)

// attr(x, '.CLASS') or CLASS(x)
#define  GET_xtsCLASS(x)                getAttrib(x, xts_ClassSymbol)
#define  SET_xtsCLASS(x,value)          setAttrib(x, xts_ClassSymbol, value)

/*
IMPORTS FROM zoo
*/
SEXP(*zoo_lag)(SEXP,SEXP,SEXP);
SEXP(*zoo_coredata)(SEXP,SEXP);

/*
FUNCTIONS
*/
SEXP do_xtsAttributes(SEXP x);              // xtsAttributes i.e. user-added attributes
SEXP do_xtsCoreAttributes(SEXP x);          /* xtsCoreAttributes xts-specific attributes
                                               CLASS, .indexFORMAT, .indexCLASS & class */
SEXP coredata(SEXP x, SEXP copyAttr);
SEXP coredata_xts(SEXP x);
SEXP add_class(SEXP x, SEXP klass);
SEXP lagXts(SEXP x, SEXP k, SEXP pad);
SEXP do_is_ordered(SEXP x, SEXP increasing, SEXP strictly);
SEXP mergeXts(SEXP args);
SEXP do_rbind_xts(SEXP x, SEXP y, SEXP dup);
SEXP rbindXts(SEXP args);
SEXP do_subset_xts(SEXP x, SEXP sr, SEXP sc, SEXP drop);
SEXP number_of_cols(SEXP args);
SEXP naCheck(SEXP x, SEXP check);

SEXP make_index_unique(SEXP x, SEXP eps);
SEXP make_unique(SEXP X, SEXP eps);
SEXP endpoints(SEXP x, SEXP on, SEXP addlast);
SEXP do_merge_xts(SEXP x, SEXP y, SEXP all, SEXP fill, SEXP retclass, SEXP colnames, 
                  SEXP suffixes, SEXP retside, SEXP env, int coerce);
SEXP na_omit_xts(SEXP x);
SEXP na_locf(SEXP x, SEXP fromlast, SEXP maxgap, SEXP limit);

SEXP tryXts(SEXP x);

SEXP xts_period_min(SEXP data, SEXP index);
SEXP xts_period_max(SEXP data, SEXP index);
SEXP xts_period_sum(SEXP data, SEXP index);
SEXP xts_period_prod(SEXP data, SEXP index);

SEXP xts_set_dimnames(SEXP x, SEXP value);


void copyAttributes(SEXP x, SEXP y);    // internal only
void copy_xtsAttributes(SEXP x, SEXP y);    // internal only
void copy_xtsCoreAttributes(SEXP x, SEXP y);// internal only    

int isXts(SEXP x);                          // is.xts analogue
int firstNonNA(SEXP x);
SEXP extract_col (SEXP x, SEXP j, SEXP drop, SEXP first_, SEXP last_);
#endif /* _XTS */

#ifdef __cplusplus
}
#endif
