/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner, Carol Woodward, Slaven Peles @ LLNL
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
 * This is the header file for a generic linear solver package.
 * It defines the SUNLinearSolver structure (_generic_SUNLinearSolver)
 * which contains the following fields:
 *   - an implementation-dependent 'content' field which contains
 *     any internal data required by the solver
 *   - an 'ops' filed which contains a structure listing operations
 *     acting on/by such solvers
 *
 * We consider both direct linear solvers and iterative linear solvers
 * as available implementations of this package; as a result some of 
 * the routines are applicable only to one type of linear solver (as 
 * noted in the comments below).
 *
 * For most of the iterative linear solvers, instead of solving the 
 * linear system A x = b directly, we apply the underlying iterative 
 * algorithm to the transformed system 
 *
 *        Abar xbar = bbar
 *
 * where
 *
 *        Abar = S1 (P1-inverse) A (P2-inverse) (S2-inverse),
 *        bbar = S1 (P1-inverse) b,
 *        xbar = S2 P2 x,
 * 
 * and where 
 *
 *        P1 = left preconditioner
 *        P2 = right preconditioner
 *        S1 = diagonal matrix of scale factors for P1-inverse b
 *        S2 = diagonal matrix of scale factors for P2 x.
 * 
 * The stopping tolerance on iterative linear solvers is on the  
 * 2-norm of the scaled preconditioned residual:
 *      || bbar - Abar xbar ||_2  <  tol.
 * 
 * We note that the Preconditioned Conjugate Gradient (PCG) solver 
 * considers S1=S2 and P1=P2 (where each is approximately A^{1/2}), and 
 * both scaling and preconditioning are applied symmetrically and 
 * simultaneously, so that the user-supplied S and P are in fact 
 * S1^2 and P1^2, i.e. P is approximately A and S is the corresponding 
 * diagonal matrix of scale factors for P.  As such, in the PCG solver the 
 * second scaling vector and the left/right "type" of P are ignored.
 *
 * -----------------------------------------------------------------
 *
 * Part I of this file contains enumeration constants for all 
 * SUNDIALS-defined linear solver types, as well as a generic type for 
 * user-supplied linear solver types.
 *
 * Part II of this file contains type declarations for the
 * _generic_SUNLinearSolver and _generic_SUNLinearSolver_Ops structures, 
 * as well as references to pointers to such structures 
 * (SUNLinearSolver).
 *
 * Part III of this file contains the prototypes for the linear solver
 * functions which operate on/by SUNLinearSolver objects.
 *
 * At a minimum, a particular implementation of a SUNLinearSolver must
 * do the following:
 *  - specify the 'content' field of SUNLinearSolver,
 *  - implement the operations on/by those SUNLinearSolver,
 *  - provide a constructor routine for new SUNLinearSolver objects
 *
 * Additionally, a SUNLinearSolver implementation may provide the 
 * following:
 *  - "Set" routines to control solver-specific parameters/options
 *  - "Get" routines to access solver-specific performance metrics
 *
 * Part IV of this file contains error codes that may be returned by
 * SUNLinearSolver objects.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINEARSOLVER_H
#define _SUNLINEARSOLVER_H
 
#include <sundials/sundials_types.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
 
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
 
  
/*
 * -----------------------------------------------------------------
 * I. Implemented SUNLinearSolver types: 
 *
 * These type names may be modified, but a at a minimum a client 
 * nonlinear solver and/or time integrator will want to know whether 
 * matrix/factorization information can be reused (hence the DIRECT 
 * and ITERATIVE types).
 * -----------------------------------------------------------------
 */
 
typedef enum {
  SUNLINEARSOLVER_DIRECT,
  SUNLINEARSOLVER_ITERATIVE,
  SUNLINEARSOLVER_CUSTOM
} SUNLinearSolver_Type;

  
/* 
 * -----------------------------------------------------------------
 * II. Generic definition of SUNLinearSolver 
 * -----------------------------------------------------------------
 */
 
/* Forward reference for pointer to SUNLinearSolver_Ops object */
typedef struct _generic_SUNLinearSolver_Ops *SUNLinearSolver_Ops;
 
/* Forward reference for pointer to SUNLinearSolver object */
typedef struct _generic_SUNLinearSolver *SUNLinearSolver;
 
/* Structure containing function pointers to linear solver operations */  
struct _generic_SUNLinearSolver_Ops {
  SUNLinearSolver_Type (*gettype)(SUNLinearSolver);
  int                  (*setatimes)(SUNLinearSolver, void*, ATimesFn);
  int                  (*setpreconditioner)(SUNLinearSolver, void*, 
                                            PSetupFn, PSolveFn);
  int                  (*setscalingvectors)(SUNLinearSolver,
                                            N_Vector, N_Vector);
  int                  (*initialize)(SUNLinearSolver);
  int                  (*setup)(SUNLinearSolver, SUNMatrix);
  int                  (*solve)(SUNLinearSolver, SUNMatrix, N_Vector, 
                                N_Vector, realtype);
  int                  (*numiters)(SUNLinearSolver);
  realtype             (*resnorm)(SUNLinearSolver);
  long int             (*lastflag)(SUNLinearSolver);
  int                  (*space)(SUNLinearSolver, long int*, long int*);
  N_Vector             (*resid)(SUNLinearSolver);
  int                  (*free)(SUNLinearSolver);
};
 
/* A linear solver is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of linear solver
   operations corresponding to that implementation. */
struct _generic_SUNLinearSolver {
  void *content;
  struct _generic_SUNLinearSolver_Ops *ops;
};

  
/*
 * -----------------------------------------------------------------
 * III. Functions exported by SUNLinearSolver module
 * 
 * SUNLinSolGetType
 *   Returns an identifier for the linear solver type from 
 *   enumeration SUNLinearSolver_Type.
 *
 * SUNLinSolSetATimes (iterative methods only)
 *   Sets the function pointer for ATimes inside of an iterative 
 *   linear solver object.  This function should only be called by a 
 *   main integrator, who will either provide this via 
 *   difference-quotients and vector operations, or by translating 
 *   between the generic ATimes call and the integrator-specific, 
 *   user-supplied routines.  This should return zero for a 
 *   successful call, and a negative value for a failure.  Ideally, 
 *   this should return one of the generic SUNLS_* error codes 
 *   listed at the bottom of this file.
 *
 * SUNLinSolSetPreconditioner (iterative methods only)
 *   Sets function pointers for PSetup and PSolve routines inside 
 *   of iterative linear solver objects.  This function should only 
 *   be called by a main integrator, who will provide translation 
 *   between the generic PSetup and PSolve calls and the integrator-
 *   specific user-supplied routines.  This should return
 *   zero for a successful call, and a negative value for a failure.
 *   Ideally, this should return one of the generic SUNLS_* error 
 *   codes listed at the bottom of this file.
 *
 * SUNLinSolSetScalingVectors (iterative methods only)
 *   Sets pointers to left/right scaling vectors for the linear
 *   system solve.  Here, s1 is an N_Vector of positive scale factors
 *   for P1-inv b, where P1 is the left preconditioner (the vector is 
 *   not tested for positivity).  Pass NULL if no scaling on P1-inv b
 *   is required.  Similarly, s2 is an N_Vector of positive scale 
 *   factors for P2 x, where P2 is the right preconditioner (again 
 *   not tested for positivity). Pass NULL if no scaling on P2 x is 
 *   required.  This should return zero for a successful call, and a 
 *   negative value for a failure.  Ideally, this should return one of
 *   the generic SUNLS_* error codes listed at the bottom of this file.
 *
 * SUNLinSolInitialize
 *   Performs linear solver initialization (assumes that all 
 *   solver-specific options have been set).  This should return
 *   zero for a successful call, and a negative value for a failure.  
 *   Ideally, this should return one of the generic SUNLS_* error 
 *   codes listed at the bottom of this file.
 *
 * SUNLinSolSetup
 *   Performs any linear solver setup needed, based on an updated
 *   system matrix A.  This may be called frequently (e.g. with a
 *   full Newton method) or infrequently (for a modified Newton 
 *   method), based on the type of integrator and/or nonlinear 
 *   solver requesting the solves.  This should return
 *   zero for a successful call, a positive value for a recoverable 
 *   failure and a negative value for an unrecoverable failure.  
 *   Ideally, this should return one of the generic SUNLS_* error 
 *   codes listed at the bottom of this file.
 *
 * SUNLinSolSolve
 *   Solves a linear system A*x = b.  If the solver is scaled, it 
 *   uses the supplied scaling vectors.  If the solver is iterative, 
 *   it attempts to solve to the specified tolerance (weighted 
 *   2-norm).  If the solver is direct it ignores the input tolerance 
 *   and scaling vectors, and if the solver does not support scaling 
 *   then it should just use a 2-norm.  This should return zero for a
 *   successful call, a positive value for a recoverable failure and 
 *   a negative value for an unrecoverable failure.  Ideally, this 
 *   should return one of the generic SUNLS_* error codes listed at 
 *   the bottom of this file.
 *
 * SUNLinSolNumIters (iterative methods only)
 *   Returns the number of linear iterations performed in the last 
 *   'Solve' call.  
 *
 * SUNLinSolResNorm (iterative methods only)
 *   Returns the final residual norm from the last 'Solve' call.
 *
 * SUNLinSolResid (iterative methods only)
 *   If an iterative method computes the preconditioned initial 
 *   residual and returns with a successful solve without performing 
 *   any iterations (i.e. either the initial guess or the 
 *   preconditioner is sufficiently accurate), then this 
 *   function may be called.  It should return the N_Vector 
 *   containing the preconditioned initial residual.
 *
 * SUNLinSolLastFlag (optional)
 *   Returns the last error flag encountered within the linear solver,
 *   allowing the user to investigate linear solver issues after 
 *   failed solves.
 *
 * SUNLinSolSpace (optional)
 *   Returns the integer and real workspace sizes for the linear
 *   solver.
 *
 * SUNLinSolFree
 *   Frees memory allocated by the linear solver.  This should return
 *   zero for a successful call, and a negative value for a failure.
 *
 * ---------------------------------------------------------------
 *
 * The following table lists the linear solver functions used by 
 * different modules in SUNDIALS.  The symbols in the table have 
 * The following meaning:
 * M - called by the modified Newton nonlinear solver
 * I - called by the inexact Newton nonlinear solver
 * P - called by the Picard nonlinear solver
 * S - called by the integrator directly (and non-identity mass matrix)
 *
 *  LinearSolver
 *  Functions           CVODE(S)  ARKode     IDA(S)      KINSOL
 *  ------------------------------------------------------------
 *  GetType             M I P#    M I P#     M I P#      M I P
 *  SetATimes           I         I S+       I           I
 *  SetPreconditioner   I*        I* S*      I*          I*
 *  SetScalingVectors   I         I S+       I           I 
 *  Initialize          M I P#    M I P# S   M I P#      M I P
 *  Setup               M I P#    M I P# S   M I P#      M I P      
 *  Solve               M I P#    M I P# S   M I P#      M I P
 *  NumIters            I         I S+       I           I
 *  ResNorm             I         I S+       I           I
 *  Resid                                    I
 *  LastFlag^
 *  Space               M I P#    M I P# S   M I P#      M I P
 *  Free                M I P#    M I P# S   M I P#      M I P
 *  ------------------------------------------------------------
 * Notes: * -- only if user calls integrator-specific 
 *             preconditioner "set" routine
 *        + -- only called when using a non-identity mass matrix
 *             with an iterative linear solver
 *        # -- planned (currently only available for KINSOL)         
 *        ^ -- available for users to diagnose solver failures
 * ---------------------------------------------------------------
 */
  
SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType(SUNLinearSolver S);

SUNDIALS_EXPORT int SUNLinSolSetATimes(SUNLinearSolver S, void* A_data,
                                       ATimesFn ATimes);
  
SUNDIALS_EXPORT int SUNLinSolSetPreconditioner(SUNLinearSolver S, void* P_data,
                                               PSetupFn Pset, PSolveFn Psol);
  
SUNDIALS_EXPORT int SUNLinSolSetScalingVectors(SUNLinearSolver S, N_Vector s1,
                                               N_Vector s2);
  
SUNDIALS_EXPORT int SUNLinSolInitialize(SUNLinearSolver S);
  
SUNDIALS_EXPORT int SUNLinSolSetup(SUNLinearSolver S, SUNMatrix A);
  
SUNDIALS_EXPORT int SUNLinSolSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                   N_Vector b, realtype tol);
  
SUNDIALS_EXPORT int SUNLinSolNumIters(SUNLinearSolver S);
  
SUNDIALS_EXPORT realtype SUNLinSolResNorm(SUNLinearSolver S);
  
SUNDIALS_EXPORT N_Vector SUNLinSolResid(SUNLinearSolver S);
  
SUNDIALS_EXPORT long int SUNLinSolLastFlag(SUNLinearSolver S);
  
SUNDIALS_EXPORT int SUNLinSolSpace(SUNLinearSolver S, long int *lenrwLS,
                                   long int *leniwLS);
  
SUNDIALS_EXPORT int SUNLinSolFree(SUNLinearSolver S);


/*
 * -----------------------------------------------------------------
 * IV. SUNLinearSolver error codes
 * ---------------------------------------------------------------
 */

#define SUNLS_SUCCESS             0  /* successful/converged          */

#define SUNLS_MEM_NULL           -1  /* mem argument is NULL          */
#define SUNLS_ILL_INPUT          -2  /* illegal function input        */
#define SUNLS_MEM_FAIL           -3  /* failed memory access          */
#define SUNLS_ATIMES_FAIL_UNREC  -4  /* atimes unrecoverable failure  */
#define SUNLS_PSET_FAIL_UNREC    -5  /* pset unrecoverable failure    */
#define SUNLS_PSOLVE_FAIL_UNREC  -6  /* psolve unrecoverable failure  */
#define SUNLS_PACKAGE_FAIL_UNREC -7  /* external package unrec. fail  */
#define SUNLS_GS_FAIL            -8  /* Gram-Schmidt failure          */        
#define SUNLS_QRSOL_FAIL         -9  /* QRsol found singular R        */

#define SUNLS_RES_REDUCED         1  /* nonconv. solve, resid reduced */
#define SUNLS_CONV_FAIL           2  /* nonconvergent solve           */
#define SUNLS_ATIMES_FAIL_REC     3  /* atimes failed recoverably     */
#define SUNLS_PSET_FAIL_REC       4  /* pset failed recoverably       */
#define SUNLS_PSOLVE_FAIL_REC     5  /* psolve failed recoverably     */
#define SUNLS_PACKAGE_FAIL_REC    6  /* external package recov. fail  */
#define SUNLS_QRFACT_FAIL         7  /* QRfact found singular matrix  */
#define SUNLS_LUFACT_FAIL         8  /* LUfact found singular matrix  */

#ifdef __cplusplus
}
#endif
#endif
