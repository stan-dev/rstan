/*
 * -----------------------------------------------------------------
 * $Revision: 4895 $
 * $Date: 2016-09-05 16:57:24 -0700 (Mon, 05 Sep 2016) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 *------------------------------------------------------------------
 * SUNDIALS configuration header file
 *------------------------------------------------------------------
 */

/* Define SUNDIALS version number */
#define SUNDIALS_VERSION_MAJOR 3
#define SUNDIALS_VERSION_MINOR 1
#define SUNDIALS_VERSION_PATCH 0

#define SUNDIALS_VERSION "3.1.0"
#define SUNDIALS_VERSION_LABEL ""


/* FCMIX: Define Fortran name-mangling macro for C identifiers.
 * Depending on the inferred scheme, one of the following six
 * macros will be defined:
 *     #define SUNDIALS_F77_FUNC(name,NAME) name
 *     #define SUNDIALS_F77_FUNC(name,NAME) name ## _
 *     #define SUNDIALS_F77_FUNC(name,NAME) name ## __
 *     #define SUNDIALS_F77_FUNC(name,NAME) NAME
 *     #define SUNDIALS_F77_FUNC(name,NAME) NAME ## _
 *     #define SUNDIALS_F77_FUNC(name,NAME) NAME ## __
 */


/* FCMIX: Define Fortran name-mangling macro for C identifiers
 *        which contain underscores.
 */


/* Define precision of SUNDIALS data type 'realtype' 
 * Depending on the precision level, one of the following 
 * three macros will be defined:
 *     #define SUNDIALS_SINGLE_PRECISION 1
 *     #define SUNDIALS_DOUBLE_PRECISION 1
 *     #define SUNDIALS_EXTENDED_PRECISION 1
 */
#define SUNDIALS_DOUBLE_PRECISION 1

/* Use generic math functions 
 * If it was decided that generic math functions can be used, then
 *     #define SUNDIALS_USE_GENERIC_MATH
 */
#define SUNDIALS_USE_GENERIC_MATH

/* Use POSIX timers if available.
 *     #define SUNDIALS_HAVE_POSIX_TIMERS
 */
/* #undef SUNDIALS_HAVE_POSIX_TIMERS */

/* Blas/Lapack available
 * If working libraries for Blas/lapack support were found, then
 *     #define SUNDIALS_BLAS_LAPACK
 */
/* #undef SUNDIALS_BLAS_LAPACK */

/* SUPERLUMT available
 * If working libraries for SUPERLUMT support were found, then
 *     #define SUNDIALS_SUPERLUMT 
 */
/* #undef SUNDIALS_SUPERLUMT */
/* #undef SUNDIALS_SUPERLUMT_THREAD_TYPE */

/* KLU available
 * If working libraries for KLU support were found, then
 *     #define SUNDIALS_KLU 
 */
/* #undef SUNDIALS_KLU */

/* FNVECTOR: Allow user to specify different MPI communicator
 * If it was found that the MPI implementation supports MPI_Comm_f2c, then
 *      #define SUNDIALS_MPI_COMM_F2C 1
 * otherwise
 *      #define SUNDIALS_MPI_COMM_F2C 0
 */


/* Mark SUNDIALS API functions for export/import
 * When building shared SUNDIALS libraries under Windows, use
 *      #define SUNDIALS_EXPORT __declspec(dllexport)
 * When linking to shared SUNDIALS libraries under Windows, use
 *      #define SUNDIALS_EXPORT __declspec(dllimport)
 * In all other cases (other platforms or static libraries under
 * Windows), the SUNDIALS_EXPORT macro is empty
 */
#define SUNDIALS_EXPORT
