/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on sundials_klu_impl.h and arkode_klu.h/cvode_klu.h/... 
 *     code, written by Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the KLU implementation of the 
 * SUNLINSOL module.
 * 
 * Part I contains declarations specific to the KLU implementation
 * of the supplied SUNLINSOL module.
 * 
 * Part II contains the prototype for the constructor 
 * SUNKLU as well as implementation-specific prototypes 
 * for various useful solver operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can 
 *     be found in the header file sundials_linearsolver.h.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_KLU_H
#define _SUNLINSOL_KLU_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "klu.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default KLU solver parameters */
#define SUNKLU_ORDERING_DEFAULT  1    /* COLAMD */

/* Interfaces to match 'sunindextype' with the correct KLU types/functions */
#if defined(SUNDIALS_INT64_T)
#define sun_klu_symbolic      klu_l_symbolic
#define sun_klu_numeric       klu_l_numeric
#define sun_klu_common        klu_l_common
#define sun_klu_analyze       klu_l_analyze
#define sun_klu_factor        klu_l_factor
#define sun_klu_refactor      klu_l_refactor
#define sun_klu_rcond         klu_l_rcond
#define sun_klu_condest       klu_l_condest
#define sun_klu_defaults      klu_l_defaults
#define sun_klu_free_symbolic klu_l_free_symbolic
#define sun_klu_free_numeric  klu_l_free_numeric
#elif defined(SUNDIALS_INT32_T)
#define sun_klu_symbolic      klu_symbolic
#define sun_klu_numeric       klu_numeric
#define sun_klu_common        klu_common
#define sun_klu_analyze       klu_analyze
#define sun_klu_factor        klu_factor
#define sun_klu_refactor      klu_refactor
#define sun_klu_rcond         klu_rcond
#define sun_klu_condest       klu_condest
#define sun_klu_defaults      klu_defaults
#define sun_klu_free_symbolic klu_free_symbolic
#define sun_klu_free_numeric  klu_free_numeric
#else  /* incompatible sunindextype for KLU */
#error  Incompatible sunindextype for KLU
#endif

#if defined(SUNDIALS_DOUBLE_PRECISION)
#else
#error  Incompatible realtype for KLU
#endif

/*
 * -----------------------------------------------------------------
 * PART I: KLU implementation of SUNLinearSolver
 *
 * The KLU implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     last_flag -- last error return flag from internal setup/solve
 *     first_factorize -- flag indicating whether the factorization 
 *       has ever been performed
 *     Symbolic -- KLU storage structure for symbolic 
 *       factorization components
 *     Numeric -- KLU storage structure for numeric factorization
 *        components
 *     Common -- storage structure for common KLU solver 
 *        components
 *     klu_solver -- ptr to KLU function to handle CSR/CSC
 * -----------------------------------------------------------------
 */
  
struct _SUNLinearSolverContent_KLU {
  long int         last_flag;
  int              first_factorize;
  sun_klu_symbolic *symbolic;
  sun_klu_numeric  *numeric;
  sun_klu_common   common;
  sunindextype     (*klu_solver)(sun_klu_symbolic*, sun_klu_numeric*,
                                 sunindextype, sunindextype,
                                 double*, sun_klu_common*);
};

typedef struct _SUNLinearSolverContent_KLU *SUNLinearSolverContent_KLU;

  
/*
 * -----------------------------------------------------------------
 * PART II: functions exported by sunlinsol_klu
 * 
 * CONSTRUCTOR:
 *    SUNKLU creates and allocates memory for a KLU sparse-direct 
 *      linear solver
 *
 * OTHER:
 *    SUNKLUReInit reinitializes memory and flags for a new 
 *      factorization (symbolic and numeric) to be conducted at the 
 *      next solver setup call.  This routine is useful in the 
 *      cases where the number of nonzeroes has changed or if the 
 *      structure of the linear system has changed which would 
 *      require a new symbolic (and numeric factorization).
 *
 *      The reinit_type argument governs the level of 
 *      reinitialization:
 *
 *      reinit_type = 1: The Jacobian matrix will be destroyed and 
 *                       a new one will be allocated based on the 
 *                       nnz value passed to this call. New 
 *                       symbolic and numeric factorizations will 
 *                       be completed at the next solver setup.
 *
 *      reinit_type = 2: Only symbolic and numeric factorizations 
 *                       will be completed.  It is assumed that the 
 *                       Jacobian size has not exceeded the size of 
 *                       nnz given in the sparse matrix provided to
 *                       the original constructor routine (or the 
 *                       previous SUNKLUReInit call) 
 *
 *      This routine assumes no other changes to solver use are 
 *      necessary.
 *
 *    SUNKLUSetOrdering sets the ordering used by KLU for reducing 
 *      fill in the linear solve.  Options for ordering_choice are: 
 *          0 for AMD, 
 *          1 for COLAMD, and 
 *          2 for the natural ordering.
 *      The default is 1 for COLAMD.
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNKLU(N_Vector y, SUNMatrix A);

SUNDIALS_EXPORT int SUNKLUReInit(SUNLinearSolver S, SUNMatrix A,
                                 sunindextype nnz, int reinit_type);

SUNDIALS_EXPORT int SUNKLUSetOrdering(SUNLinearSolver S,
                                      int ordering_choice);

/*
 * -----------------------------------------------------------------
 * KLU implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_KLU(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_KLU(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_KLU(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_KLU(SUNLinearSolver S, SUNMatrix A,
                                       N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_KLU(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_KLU(SUNLinearSolver S,
                                       long int *lenrwLS,
                                       long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_KLU(SUNLinearSolver S);
  

#ifdef __cplusplus
}
#endif

#endif
