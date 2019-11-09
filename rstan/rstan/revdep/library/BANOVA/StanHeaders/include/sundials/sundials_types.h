/* -----------------------------------------------------------------
 * Programmer(s): Scott Cohen, Alan Hindmarsh, Radu Serban,
 *                Aaron Collier, and Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This header file exports three types: realtype, sunindextype and 
 * booleantype, as well as the constants SUNTRUE and SUNFALSE.
 *
 * Users should include the header file sundials_types.h in every
 * program file and use the exported name realtype instead of
 * float, double or long double.
 *
 * The constants SUNDIALS_SINGLE_PRECISION, SUNDIALS_DOUBLE_PRECISION
 * and SUNDIALS_LONG_DOUBLE_PRECISION indicate the underlying data
 * type of realtype. 
 * 
 * The legal types for realtype are float, double and long double.
 *
 * The constants SUNDIALS_INT64_T and SUNDIALS_INT32_T indicate 
 * the underlying data type of sunindextype -- the integer data type
 * used for vector and matrix indices. 
 * 
 * Data types are set at the configuration stage.
 *
 * The macro RCONST gives the user a convenient way to define
 * real-valued literal constants. To use the constant 1.0, for example,
 * the user should write the following:
 *
 *   #define ONE RCONST(1.0)
 *
 * If realtype is defined as a double, then RCONST(1.0) expands
 * to 1.0. If realtype is defined as a float, then RCONST(1.0)
 * expands to 1.0F. If realtype is defined as a long double,
 * then RCONST(1.0) expands to 1.0L. There is never a need to
 * explicitly cast 1.0 to (realtype). The macro can be used for 
 * literal constants only. It cannot be used for expressions. 
 * -----------------------------------------------------------------*/
  
#ifndef _SUNDIALSTYPES_H
#define _SUNDIALSTYPES_H

#ifndef _SUNDIALS_CONFIG_H
#define _SUNDIALS_CONFIG_H
#include <sundials/sundials_config.h>
#endif

#include <float.h>
#include <stdint.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 *------------------------------------------------------------------
 * Type realtype
 * Macro RCONST
 * Constants BIG_REAL, SMALL_REAL, and UNIT_ROUNDOFF
 *------------------------------------------------------------------
 */

#if defined(SUNDIALS_SINGLE_PRECISION)

typedef float realtype;
# define RCONST(x) x##F
# define BIG_REAL FLT_MAX
# define SMALL_REAL FLT_MIN
# define UNIT_ROUNDOFF FLT_EPSILON

#elif defined(SUNDIALS_DOUBLE_PRECISION)

typedef double realtype;
# define RCONST(x) x
# define BIG_REAL DBL_MAX
# define SMALL_REAL DBL_MIN
# define UNIT_ROUNDOFF DBL_EPSILON

#elif defined(SUNDIALS_EXTENDED_PRECISION)

typedef long double realtype;
# define RCONST(x) x##L
# define BIG_REAL LDBL_MAX
# define SMALL_REAL LDBL_MIN
# define UNIT_ROUNDOFF LDBL_EPSILON

#endif


/*
 *------------------------------------------------------------------
 * Type : sunindextype
 *------------------------------------------------------------------
 * Defines integer type to be used for vector and matrix indices.
 * User can build sundials to use 32- or 64-bit signed integers.
 * If compiler does not support portable data types, 'int' is used
 * in place of 32- and 'long int' is used in place of 64-bit 
 * integer.
 *------------------------------------------------------------------
 */

#define SUNDIALS_INT64_T

#if defined(SUNDIALS_INT64_T)

#if __STDC_VERSION__ >= 199901L
typedef int64_t sunindextype;
#else
typedef long int sunindextype;
#endif

#elif defined(SUNDIALS_INT32_T)

#if __STDC_VERSION__ >= 199901L
typedef int32_t sunindextype;
#else
typedef int sunindextype;
#endif

#endif


/*
 *------------------------------------------------------------------
 * Type : booleantype
 *------------------------------------------------------------------
 * Constants : SUNFALSE and SUNTRUE
 *------------------------------------------------------------------
 * ANSI C does not have a built-in boolean data type. Below is the
 * definition for a new type called booleantype. The advantage of
 * using the name booleantype (instead of int) is an increase in
 * code readability. It also allows the programmer to make a
 * distinction between int and boolean data. Variables of type
 * booleantype are intended to have only the two values SUNFALSE and
 * SUNTRUE which are defined below to be equal to 0 and 1,
 * respectively.
 *------------------------------------------------------------------
 */

#ifndef booleantype
#define booleantype int
#endif

#ifndef SUNFALSE
#define SUNFALSE 0
#endif

#ifndef SUNTRUE
#define SUNTRUE 1
#endif


#ifdef __cplusplus
}
#endif

#endif  /* _SUNDIALSTYPES_H */


