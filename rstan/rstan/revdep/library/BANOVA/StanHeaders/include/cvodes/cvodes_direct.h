/* ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Radu Serban @ LLNL
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
 * Common header file for the direct linear solver interface in 
 * CVODES.
 *
 * Part I contains type definitions and function prototypes for 
 * using a CVDLS linear solver on forward problems (IVP 
 * integration and/or FSA)
 *
 * Part II contains type definitions and function prototypes for 
 * using a CVDLS linear solver on adjoint (backward) problems
 * -----------------------------------------------------------------*/

#ifndef _CVSDLS_H
#define _CVSDLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  CVSDLS Constants
  =================================================================*/

  
/*-----------------------------------------------------------------
  return values 
  -----------------------------------------------------------------*/

#define CVDLS_SUCCESS           0
#define CVDLS_MEM_NULL         -1
#define CVDLS_LMEM_NULL        -2
#define CVDLS_ILL_INPUT        -3
#define CVDLS_MEM_FAIL         -4

/* Additional last_flag values */

#define CVDLS_JACFUNC_UNRECVR  -5
#define CVDLS_JACFUNC_RECVR    -6
#define CVDLS_SUNMAT_FAIL      -7

/* Return values for the adjoint module */

#define CVDLS_NO_ADJ           -101
#define CVDLS_LMEMB_NULL       -102


/*=================================================================
  PART I:  Forward Problems
  =================================================================*/

/*-----------------------------------------------------------------
  FUNCTION TYPES
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  Type: CVDlsJacFn
  -----------------------------------------------------------------
 
  A Jacobian approximation function Jac must be of type CVDlsJacFn. 
  Its parameters are:
 
  Jac is the SUNMatrix that will be loaded by a CVDlsJacFn with an 
      approximation to the Jacobian matrix J = (df_i/dy_j) at 
      the point (t,y). 
 
  t   is the current value of the independent variable.
 
  y   is the current value of the dependent variable vector,
      namely the predicted value of y(t).
 
  fy  is the vector f(t,y).
 
  user_data is a pointer to user data - the same as the user_data
      parameter passed to CVodeSetUserdata.
 
  tmp1, tmp2, and tmp3 are pointers to memory allocated for
  vectors of length N which can be used by a CVDlsJacFn
  as temporary storage or work space.
 
  A CVDlsJacFn should return 0 if successful, a positive 
  value if a recoverable error occurred, and a negative value if 
  an unrecoverable error occurred.
 
  -----------------------------------------------------------------

  NOTE: See the relevant SUNMatrix implementation header files and
      documentation for mechanisms to inquire about matrix 
      dimensions, and for efficient ways to set matrix entries.
                                                                 
  NOTE: If the user's Jacobian routine needs other quantities,   
      they are accessible as follows: hcur (the current stepsize)
      and ewt (the error weight vector) are accessible through   
      CVodeGetCurrentStep and CVodeGetErrWeights, respectively 
      (see cvode.h). The unit roundoff is available as 
      UNIT_ROUNDOFF defined in sundials_types.h.
 
  -----------------------------------------------------------------*/
typedef int (*CVDlsJacFn)(realtype t, N_Vector y, N_Vector fy, 
                          SUNMatrix Jac, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  

/*-----------------------------------------------------------------
  EXPORTED FUNCTIONS 
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  Required inputs to the CVSDLS linear solver interface
  -----------------------------------------------------------------
 
  CVDlsSetLinearSolver specifies the direct SUNLinearSolver object
  that should be used.  This is required if CVodes is solving a 
  problem using the Newton nonlinear solver (not the function 
  iteration).
 
  The return value is one of:
     CVDLS_SUCCESS   if successful
     CVDLS_MEM_NULL  if the CVODE memory was NULL
     CVDLS_ILL_INPUT if the arguments are incompatible
  -----------------------------------------------------------------*/
SUNDIALS_EXPORT int CVDlsSetLinearSolver(void *cvode_mem, 
                                         SUNLinearSolver LS,
                                         SUNMatrix A);

  
/*-----------------------------------------------------------------
  Optional inputs to the CVSDLS linear solver
  -----------------------------------------------------------------
 
  CVDlsSetJacFn specifies the dense/band/sparse Jacobian
  approximation routine to be used for a direct linear solver.
  By default, a difference quotient approximation is used for 
  dense/band; no default exists for sparse (so this must be 
  user-supplied).
 
  The return value is one of:
     CVDLS_SUCCESS   if successful
     CVDLS_MEM_NULL  if the CVODE memory was NULL
     CVDLS_LMEM_NULL if the linear solver memory was NULL
  -----------------------------------------------------------------*/
SUNDIALS_EXPORT int CVDlsSetJacFn(void *cvode_mem, CVDlsJacFn jac);


/*-----------------------------------------------------------------
  Optional outputs from the CVSSDLS linear solver interface
  -----------------------------------------------------------------
 
  CVDlsGetWorkSpace   returns the real and integer workspace used
                      by the direct linear solver.
  CVDlsGetNumJacEvals returns the number of calls made to the
                      Jacobian evaluation routine jac.
  CVDlsGetNumRhsEvals returns the number of calls to the user
                      f routine due to finite difference Jacobian
                      evaluation.
  CVDlsGetLastFlag    returns the last error flag set by any of
                      the CVSDIRECT interface functions.
 
  The return value of CVDlsGet* is one of:
     CVDLS_SUCCESS   if successful
     CVDLS_MEM_NULL  if the CVODES memory was NULL
     CVDLS_LMEM_NULL if the linear solver memory was NULL
  -----------------------------------------------------------------*/

SUNDIALS_EXPORT int CVDlsGetWorkSpace(void *cvode_mem,
                                      long int *lenrwLS,
                                      long int *leniwLS);
SUNDIALS_EXPORT int CVDlsGetNumJacEvals(void *cvode_mem,
                                        long int *njevals);
SUNDIALS_EXPORT int CVDlsGetNumRhsEvals(void *cvode_mem,
                                        long int *nfevalsLS);
SUNDIALS_EXPORT int CVDlsGetLastFlag(void *cvode_mem,
                                     long int *flag);

/*-----------------------------------------------------------------
  The following function returns the name of the constant 
  associated with a CVSDLS return flag
  -----------------------------------------------------------------*/
SUNDIALS_EXPORT char *CVDlsGetReturnFlagName(long int flag);

/*=================================================================
  PART II:  Backward Problems
  =================================================================*/

/*-----------------------------------------------------------------
  FUNCTION TYPES
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  Type: CVDlsJacFnB
  -----------------------------------------------------------------
  A Jacobian approximation function jacB for the adjoint
  (backward) problem must have the prototype given below. 
  -----------------------------------------------------------------*/
typedef int (*CVDlsJacFnB)(realtype t, N_Vector y, N_Vector yB,
                           N_Vector fyB, SUNMatrix JB,
                           void *user_dataB, N_Vector tmp1B,
                           N_Vector tmp2B, N_Vector tmp3B);

/*-----------------------------------------------------------------
  Type: CVDlsJacFnBS
  -----------------------------------------------------------------
  A Jacobian approximation function jacBS for the adjoint
  (backward) problem, sensitivity-dependent case,  must have the
  prototype given below. 
  -----------------------------------------------------------------*/
typedef int (*CVDlsJacFnBS)(realtype t, N_Vector y, N_Vector *yS,
                            N_Vector yB, N_Vector fyB, SUNMatrix JB,
                            void *user_dataB, N_Vector tmp1B,
                            N_Vector tmp2B, N_Vector tmp3B);


/*-----------------------------------------------------------------
  EXPORTED FUNCTIONS 
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  Required inputs for the CVSDLS linear solver interface:

  CVDlsSetLinearSolverB specifies the direct SUNLinearSolver 
  object that should be used for backward integration. The 
  'which' argument is the int returned by CVodeCreateB.

  The return value is one of:
    CVDLS_SUCCESS   if successful
    CVDLS_MEM_NULL  if the cvode memory was NULL
    CVDLS_ILL_INPUT if the arguments are incompatible
---------------------------------------------------------------*/
SUNDIALS_EXPORT int CVDlsSetLinearSolverB(void *cvode_mem,
                                          int which,
                                          SUNLinearSolver LS,
                                          SUNMatrix A);
  
/*-----------------------------------------------------------------
  Functions: CVDlsSetJacFnB and CVDlsSetJacFnBS
  -----------------------------------------------------------------
  CVDlsSetJacFnB specifies the Jacobian function to be used by a
  CVSDLS linear solver for the backward integration phase, when
  the backward problem does not depend on forward sensitivities.
  CVDlsSetJacFnBS specifies the Jacobian function when the backward
  problem does depend on sensitivities.

  The 'which' argument is the int returned by CVodeCreateB.
 
  The return value is one of:
     CVDLS_SUCCESS   if successful
     CVDLS_MEM_NULL  if the CVODE memory was NULL
     CVDLS_LMEM_NULL if the linear solver memory was NULL
  -----------------------------------------------------------------*/
SUNDIALS_EXPORT int CVDlsSetJacFnB(void *cvode_mem, int which,
                                   CVDlsJacFnB jacB);
SUNDIALS_EXPORT int CVDlsSetJacFnBS(void *cvode_mem, int which,
                                    CVDlsJacFnBS jacBS);

#ifdef __cplusplus
}
#endif

#endif
