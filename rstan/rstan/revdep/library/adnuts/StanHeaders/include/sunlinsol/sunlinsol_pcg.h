/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds, Ashley Crawford @ SMU
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
 * This is the header file for the PCG implementation of the 
 * SUNLINSOL module.  The PCG algorithm is based on the
 * Preconditioned Conjugate Gradient.
 *
 * The PCG algorithm solves a linear system A x = b where 
 * A is a symmetric, real-valued matrix, i.e. A = A' (Matlab 
 * notation for the transpose of A).  Preconditioning is allowed, 
 * and is applied in a symmetric fashion on both the right and left.  
 * Scaling is also allowed and is applied symmetrically.  We denote
 * the preconditioner and scaling matrices as follows:
 *   P = preconditioner (assumed symmetric)
 *   S = diagonal matrix of scale factors 
 * The matrices A and P are not required explicitly; only routines 
 * that provide A and P-inverse as operators are required.  The 
 * diagonal of the matrix S is held in a single N_Vector, supplied
 * by the user of this module.
 *
 * In this notation, PCG applies the underlying algorithm to
 * the equivalent transformed system
 *   Abar xbar = bbar , where
 *   Abar = S (P-inverse) A (P-inverse) S ,
 *   bbar = S (P-inverse) b , and   xbar = (S-inverse) P x .
 *
 * The scaling matrix must be chosen so that the vectors
 * (S P-inverse b) and (S-inverse P x) have dimensionless 
 * components.
 *
 * The stopping test for the PCG iterations is on the L2 norm of
 * the scaled preconditioned residual:
 *      || bbar - Abar xbar ||_2  <  delta
 *  <=>
 *      || S (P-inverse) b - S (P-inverse) A x ||_2  <  delta
 *  <=>
 *      || P-inverse b - (P-inverse) A x ||_S  <  delta
 * where || v ||_S =  sqrt(v' S' S v) with an input test constant 
 * delta.
 *
 * The usage of this PCG solver involves supplying up to three
 * routines and making a variety of calls.  The user-supplied 
 * routines are
 *    atimes (A_data, x, y) to compute y = A x, given x,
 *    psolve (P_data, y, x, lr) to solve P1 x = y or P2 x = y for 
 *           x, given y,
 *    psetup (P_data) to perform any 'setup' operations in 
 *           preparation for calling psolve.
 * The user calls are:
 *    SUNLinearSolver LS = SUNPCG(y, pretype, maxl);
 *           to create the linear solver structure,
 *    flag = SUNLinSolSetATimes(LS, A_data, atimes);
 *           to set the matrix-vector product setup/apply routines,
 *    flag = SUNLinSolSetPreconditioner(LS, P_data, psetup, psolve);
 *           to *optionally* set the preconditioner setup/apply routines,
 *    flag = SUNLinSolSetScalingVectors(LS, s, NULL);
 *           to *optionally* set the diagonal of the scaling matrix 
 *           (for PCG, only the first of the two scaling vectors is used)
 *    flag = SUNLinSolInitialize(LS);
 *           to perform internal solver memory allocations,
 *    flag = SUNLinSolSetup(LS, NULL);
 *           to call the psetup routine (if non-NULL);
 *    flag = SUNLinSolSolve(LS, NULL, x, b, w, tol);
 *           to solve the linear system to the tolerance 'tol'
 *    long int nli = SUNLinSolNumIters(LS);
 *           to *optionally* retrieve the number of linear iterations 
 *           performed by the solver,
 *    long int lastflag = SUNLinSolLastFlag(LS);
 *           to *optionally* retrieve the last internal solver error flag,
 *    realtype resnorm = SUNLinSolResNorm(LS);
 *           to *optionally* retrieve the final linear residual norm,
 *    flag = SUNLinSolFree(LS);
 *           to free the solver memory.
 * Complete details for specifying atimes, psetup and psolve 
 * and for the usage calls are given below.
 *
 * -----------------------------------------------------------------
 * 
 * Part I contains declarations specific to the PCG implementation
 * of the supplied SUNLINSOL module.
 * 
 * Part II contains the prototype for the constructor SUNPCG as well
 * as implementation-specific prototypes for various useful solver
 * operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can 
 *     be found in the header file sundials_linearsolver.h.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_PCG_H
#define _SUNLINSOL_PCG_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_pcg.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default PCG solver parameters */
#define SUNPCG_MAXL_DEFAULT    5

/*
 * -----------------------------------------------------------------
 * PART I: PCG implementation of SUNLinearSolver
 *
 * The PCG implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     maxl -- number of PCG iterations to use
 *     pretype -- flag for use of preconditioning
 *     numiters -- number of iterations from most-recent solve
 *     resnorm -- final linear residual norm from most-recent solve
 *     last_flag -- last error return flag from internal setup/solve
 *     ATimes -- function pointer to ATimes routine
 *     ATData -- pointer to structure for ATimes
 *     Psetup -- function pointer to preconditioner setup routine
 *     Psolve -- function pointer to preconditioner solve routine
 *     PData -- pointer to structure for Psetup/Psolve
 *     s -- vector (type N_Vector) which holds the diagonal of the 
 *         scaling matrix S
 *     r -- vector (type N_Vector) which holds the preconditioned 
 *         linear system residual
 *     p, z and Ap -- vectors (type N_Vector) used for workspace by
 *         the PCG algorithm
 * -----------------------------------------------------------------
 */
 
struct _SUNLinearSolverContent_PCG {
  int maxl;
  int pretype;
  int numiters;
  realtype resnorm;
  long int last_flag;

  ATimesFn ATimes;
  void* ATData;
  PSetupFn Psetup;
  PSolveFn Psolve;
  void* PData;

  N_Vector s;
  N_Vector r;
  N_Vector p;
  N_Vector z;
  N_Vector Ap;
};

typedef struct _SUNLinearSolverContent_PCG *SUNLinearSolverContent_PCG;

  
/*
 * -----------------------------------------------------------------
 * PART III: functions exported by sunlinsol_pcg
 * 
 * CONSTRUCTOR:
 *    SUNPCG creates and allocates memory for a PCG solver
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNPCG(N_Vector y, int pretype, int maxl);
SUNDIALS_EXPORT int SUNPCGSetPrecType(SUNLinearSolver S, int pretype);
SUNDIALS_EXPORT int SUNPCGSetMaxl(SUNLinearSolver S, int maxl);

/*
 * -----------------------------------------------------------------
 * PCG implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_PCG(SUNLinearSolver S, void* A_data,
                                           ATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_PCG(SUNLinearSolver S,
                                                   void* P_data,
                                                   PSetupFn Pset,
                                                   PSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_PCG(SUNLinearSolver S,
                                                   N_Vector s,
                                                   N_Vector nul);
SUNDIALS_EXPORT int SUNLinSolSetup_PCG(SUNLinearSolver S, SUNMatrix nul);
SUNDIALS_EXPORT int SUNLinSolSolve_PCG(SUNLinearSolver S, SUNMatrix nul,
                                       N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_PCG(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_PCG(SUNLinearSolver S, 
                                       long int *lenrwLS, 
                                       long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_PCG(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
