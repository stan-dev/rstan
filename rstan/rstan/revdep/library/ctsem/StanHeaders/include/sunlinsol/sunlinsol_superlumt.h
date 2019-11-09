/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on codes sundials_superlumt_impl.h and <solver>_superlumt.h
 *     written by Carol S. Woodward @ LLNL
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
 * This is the header file for the SuperLUMT implementation of the 
 * SUNLINSOL module.
 * 
 * Part I contains declarations specific to the SuperLUMT 
 * implementation of the supplied SUNLINSOL module.
 * 
 * Part II contains the prototype for the constructor 
 * SUNSuperLUMT as well as implementation-specific prototypes 
 * for various useful solver operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can 
 *     be found in the header file sundials_linearsolver.h.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_SLUMT_H
#define _SUNLINSOL_SLUMT_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h>

/* assume SuperLU_MT library was built with compatible index type */  
#if defined(SUNDIALS_INT64_T)
#define _LONGINT
#endif

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SuperLU_MT solver parameters */
#define SUNSLUMT_ORDERING_DEFAULT  3     /* COLAMD */

/* Interfaces to match 'realtype' with the correct SuperLUMT functions */
#if defined(SUNDIALS_DOUBLE_PRECISION)
#ifndef _SLUMT_H
#define _SLUMT_H
#include "slu_mt_ddefs.h"
#endif
#define xgstrs                  dgstrs
#define pxgstrf                 pdgstrf
#define pxgstrf_init            pdgstrf_init
#define xCreate_Dense_Matrix    dCreate_Dense_Matrix
#define xCreate_CompCol_Matrix  dCreate_CompCol_Matrix
#elif defined(SUNDIALS_SINGLE_PRECISION)
#ifndef _SLUMT_H
#define _SLUMT_H
#include "slu_mt_sdefs.h"
#endif
#define xgstrs                  sgstrs
#define pxgstrf                 psgstrf
#define pxgstrf_init            psgstrf_init
#define xCreate_Dense_Matrix    sCreate_Dense_Matrix
#define xCreate_CompCol_Matrix  sCreate_CompCol_Matrix
#else  /* incompatible sunindextype for SuperLUMT */
#error  Incompatible realtype for SuperLUMT
#endif

  
/*
 * -----------------------------------------------------------------
 * PART I: SuperLUMT implementation of SUNLinearSolver
 *
 * The SuperLUMT implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     last_flag -- last error return flag from internal setup/solve
 *     first_factorize -- flag indicating whether the factorization 
 *       has ever been performed
 *     A, AC, L, U, B -- SuperMatrix pointers used in solve
 *     Gstat -- GStat_t object used in solve
 *     perm_r, perm_c -- permutation arrays used in solve
 *     N -- size of the linear system
 *     num_threads -- number of OpenMP/Pthreads threads to use
 *     diag_pivot_thresh -- threshold on diagonal pivoting
 *     ordering -- flag for reordering algorithm to use
 *     options -- pointer to SuperLUMT options structure
 * -----------------------------------------------------------------
 */
  
struct _SUNLinearSolverContent_SuperLUMT {
  long int     last_flag;
  int          first_factorize;
  SuperMatrix  *A, *AC, *L, *U, *B;
  Gstat_t      *Gstat;
  sunindextype *perm_r, *perm_c;
  sunindextype N;
  int          num_threads;
  realtype     diag_pivot_thresh; 
  int          ordering;
  superlumt_options_t *options;
};

typedef struct _SUNLinearSolverContent_SuperLUMT *SUNLinearSolverContent_SuperLUMT;

  
/*
 * -----------------------------------------------------------------
 * PART II: functions exported by sunlinsol_slumt
 * 
 * CONSTRUCTOR:
 *    SUNSuperLUMT creates and allocates memory for a SuperLUMT sparse-direct 
 *      linear solver
 *
 * OTHER:
 *    SUNSuperLUMTSetOrdering sets the ordering used by SuperLUMT for reducing 
 *      fill in the linear solve.  Options for ordering_choice are: 
 *         0 for natural ordering
 *         1 for minimal degree ordering on A'*A
 *         2 for minimal degree ordering on A'+A
 *         3 for AMD ordering for unsymmetric matrices
 *      The default used in SUNDIALS is 3 for COLAMD.
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNSuperLUMT(N_Vector y, SUNMatrix A,
                                             int num_threads);

SUNDIALS_EXPORT int SUNSuperLUMTSetOrdering(SUNLinearSolver S,
                                            int ordering_choice);

/*
 * -----------------------------------------------------------------
 * SuperLUMT implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SuperLUMT(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_SuperLUMT(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_SuperLUMT(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_SuperLUMT(SUNLinearSolver S, SUNMatrix A,
                                       N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_SuperLUMT(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_SuperLUMT(SUNLinearSolver S, 
                                             long int *lenrwLS, 
                                             long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_SuperLUMT(SUNLinearSolver S);
  

#ifdef __cplusplus
}
#endif

#endif
