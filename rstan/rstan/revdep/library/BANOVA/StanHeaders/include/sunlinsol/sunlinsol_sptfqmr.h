/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on code sundials_sptfqmr.h by: Aaron Collier @ LLNL
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
 * This is the header file for the SPTFQMR implementation of the 
 * SUNLINSOL module.  The SPTFQMR algorithm is based on the
 * Scaled Preconditioned Transpose-free Quasi-Minimum Residual method.
 *
 * The SPTFQMR algorithm solves a linear system A x = b.
 * Preconditioning is allowed on the left, right, or both.
 * Scaling is allowed on both sides.
 * We denote the preconditioner and scaling matrices as follows:
 *   P1 = left preconditioner
 *   P2 = right preconditioner
 *   S1 = diagonal matrix of scale factors for P1-inverse b
 *   S2 = diagonal matrix of scale factors for P2 x
 * The matrices A, P1, and P2 are not required explicitly; only
 * routines that provide A, P1-inverse, and P2-inverse as
 * operators are required.
 *
 * In this notation, SPTFQMR applies the underlying TFQMR method to
 * the equivalent transformed system
 *   Abar xbar = bbar , where
 *   Abar = S1 (P1-inverse) A (P2-inverse) (S2-inverse) ,
 *   bbar = S1 (P1-inverse) b , and   xbar = S2 P2 x .
 *
 * The scaling matrices must be chosen so that vectors S1
 * P1-inverse b and S2 P2 x have dimensionless components.
 * If preconditioning is done on the left only (P2 = I), by a
 * matrix P, then S2 must be a scaling for x, while S1 is a
 * scaling for P-inverse b, and so may also be taken as a scaling
 * for x.  Similarly, if preconditioning is done on the right only
 * (P1 = I, P2 = P), then S1 must be a scaling for b, while S2 is
 * a scaling for P x, and may also be taken as a scaling for b.
 *
 * The stopping test for the SPTFQMR iterations is on the L2 norm of
 * the scaled preconditioned residual:
 *      || bbar - Abar xbar ||_2  <  delta
 * with an input test constant delta.
 *
 * The usage of this SPTFQMR solver involves supplying up to three 
 * routines and making a variety of calls.  The user-supplied routines are
 *    atimes (A_data, x, y) to compute y = A x, given x,
 *    psolve (P_data, y, x, lr) to solve P1 x = y or P2 x = y for 
 *           x, given y,
 *    psetup (P_data) to perform any 'setup' operations in 
 *           preparation for calling psolve.
 * The three user calls are:
 *    SUNLinearSolver LS = SUNSPTFQMR(y, pretype, maxl);
 *           to create the linear solver structure,
 *    flag = SUNLinSolSetATimes(LS, A_data, atimes);
 *           to set the matrix-vector product setup/apply routines,
 *    flag = SUNLinSolSetPreconditioner(LS, P_data, psetup, psolve);
 *           to *optionally* set the preconditioner setup/apply routines,
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
 *    flag = SUNLinSolFree(LS);
 *           to free the solver memory.
 * Complete details for specifying atimes, psetup and psolve 
 * and for the usage calls are given below.
 *
 * -----------------------------------------------------------------
 * 
 * Part I contains declarations specific to the SPTFQMR implementation
 * of the supplied SUNLINSOL module.
 * 
 * Part II contains the prototype for the constructor 
 * SUNSPTFQMR as well as implementation-specific prototypes 
 * for various useful solver operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can 
 *     be found in the header file sundials_linearsolver.h.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_SPTFQMR_H
#define _SUNLINSOL_SPTFQMR_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_sptfqmr.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Default SPTFQMR solver parameters */
#define SUNSPTFQMR_MAXL_DEFAULT    5

/*
 * -----------------------------------------------------------------
 * PART I: SPTFQMR implementation of SUNLinearSolver
 *
 * The SPTFQMR implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     maxl -- number of BiCGStab iterations to allow
 *     pretype -- flag for type of preconditioning to employ
 *     numiters -- number of iterations from most-recent solve
 *     resnorm -- final linear residual norm from most-recent solve
 *     last_flag -- last error return flag from internal setup/solve
 *     ATimes -- function pointer to ATimes routine
 *     ATData -- pointer to structure for ATimes
 *     Psetup -- function pointer to preconditioner setup routine
 *     Psolve -- function pointer to preconditioner solve routine
 *     PData -- pointer to structure for Psetup/Psolve
 *     s1, s2 -- vector pointers for supplied scaling matrices
 *     r_star -- a vector (type N_Vector) which holds the initial 
 *         scaled, preconditioned linear system residual
 *     q, d, v, p and u -- vectors (type N_Vector) used for 
 *         workspace by the SPTFQMR algorithm
 *     r -- array of vectors (type N_Vector) used for workspace 
 *         within the SPTFQMR algorithm
 *     vtemp1/vtemp2/vtemp3 -- scratch vectors (type N_Vector) used
 *         as temporary vector storage during calculations.
 * -----------------------------------------------------------------
 */
  
struct _SUNLinearSolverContent_SPTFQMR {
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

  N_Vector s1;
  N_Vector s2;
  N_Vector r_star;
  N_Vector q;
  N_Vector d;
  N_Vector v;
  N_Vector p;
  N_Vector *r;
  N_Vector u;
  N_Vector vtemp1;
  N_Vector vtemp2;
  N_Vector vtemp3;
};

typedef struct _SUNLinearSolverContent_SPTFQMR *SUNLinearSolverContent_SPTFQMR;

/*
 * -----------------------------------------------------------------
 * PART II: functions exported by sunlinsol_sptfqmr
 * 
 * CONSTRUCTOR:
 *    SUNSPTFQMR creates and allocates memory for a SPTFQMR solver
 *
 * "SET" ROUTINES:
 *    SUNSPTFQMRSSetPrecType updates the type of preconditioning to 
 *       use.  Supported values are PREC_NONE, PREC_LEFT, PREC_RIGHT 
 *       and PREC_BOTH.
 *    SUNSPTFQMRSetMaxl updates the maximum number of iterations to 
 *       allow in the solver.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNSPTFQMR(N_Vector y, int pretype, int maxl);
SUNDIALS_EXPORT int SUNSPTFQMRSetPrecType(SUNLinearSolver S, int pretype);
SUNDIALS_EXPORT int SUNSPTFQMRSetMaxl(SUNLinearSolver S, int maxl);

/*
 * -----------------------------------------------------------------
 * SPTFQMR implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetATimes_SPTFQMR(SUNLinearSolver S, void* A_data,
                                               ATimesFn ATimes);
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner_SPTFQMR(SUNLinearSolver S,
                                                       void* P_data,
                                                       PSetupFn Pset,
                                                       PSolveFn Psol);
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors_SPTFQMR(SUNLinearSolver S,
                                                       N_Vector s1,
                                                       N_Vector s2);
SUNDIALS_EXPORT int SUNLinSolSetup_SPTFQMR(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_SPTFQMR(SUNLinearSolver S, SUNMatrix A,
                                           N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT int SUNLinSolNumIters_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT realtype SUNLinSolResNorm_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT N_Vector SUNLinSolResid_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_SPTFQMR(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_SPTFQMR(SUNLinearSolver S, 
                                           long int *lenrwLS, 
                                           long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_SPTFQMR(SUNLinearSolver S);


#ifdef __cplusplus
}
#endif

#endif
