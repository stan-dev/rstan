/*
 * -----------------------------------------------------------------
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
 * Common header file for the direct linear solver interface in IDAS
 * -----------------------------------------------------------------
 */

#ifndef _IDADLS_H
#define _IDADLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  IDASDLS Constants
  =================================================================*/

/* IDASDLS return values */
#define IDADLS_SUCCESS           0
#define IDADLS_MEM_NULL         -1
#define IDADLS_LMEM_NULL        -2
#define IDADLS_ILL_INPUT        -3
#define IDADLS_MEM_FAIL         -4

/* Additional last_flag values */
#define IDADLS_JACFUNC_UNRECVR  -5
#define IDADLS_JACFUNC_RECVR    -6
#define IDADLS_SUNMAT_FAIL      -7

/* Return values for the adjoint module */
#define IDADLS_NO_ADJ           -101
#define IDADLS_LMEMB_NULL       -102

/*=================================================================
  PART I:  Forward Problems
  =================================================================*/

/*-----------------------------------------------------------------
  FUNCTION TYPES
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  Type : IDADlsJacFn
 
  A Jacobian approximation function Jac must be of type IDADlsJacFn.
  Its parameters are:                     
                                                                 
  t   is the current value of the independent variable.
                                                                 
  y   is the current value of the dependent variable vector,     
      namely the predicted value of y(t).                     
                                                                 
  yp  is the current value of the derivative vector y',          
      namely the predicted value of y'(t).                    
                                                                 
  r   is the residual vector F(tt,yy,yp).                     
                                                                 
  c_j is the scalar in the system Jacobian, proportional to 
      the inverse of the step size h.
                                                                 
  user_data is a pointer to user Jacobian data - the same as the    
      user_data parameter passed to IDASetUserData.                     
                                                                 
  Jac is the SUNMatrix to be loaded by an IDADlsJacFn routine 
      with an approximation to the system Jacobian matrix                                  
             J = dF/dy + c_j *dF/dy'
      at the given point (t,y,y'), where the ODE system is given 
      by F(t,y,y') = 0.   

      Note that Jac is NOT preset to zero!
                                                                 
  tmp1, tmp2, tmp3 are pointers to memory allocated for          
      N_Vectors which can be used by an IDADlsDenseJacFn routine 
      as temporary storage or work space.                     
                                                                 
  A IDADlsDenseJacFn should return                                
      0 if successful,                                           
      a positive int if a recoverable error occurred, or         
      a negative int if a nonrecoverable error occurred.         
  In the case of a recoverable error return, the integrator will 
  attempt to recover by reducing the stepsize (which changes c_j).
 
  NOTE: See the relevant SUNMatrix implementation header files
      and documentation for mechanisms to inquire about matrix 
      dimensions, and for efficient ways to set matrix entries.
                                                                
                                                                
  NOTE: If the user's Jacobian routine needs other quantities,   
      they are accessible as follows: hcur (the current stepsize)
      and ewt (the error weight vector) are accessible through   
      IDAGetCurrentStep and IDAGetErrWeights, respectively, but this
      requires including in user_data a pointer to the solver memory.
      The unit roundoff is available as UNIT_ROUNDOFF defined in
      sundials_types.h.
 
  -----------------------------------------------------------------*/
typedef int (*IDADlsJacFn)(realtype t, realtype c_j, N_Vector y,
                           N_Vector yp, N_Vector r, SUNMatrix Jac,
                           void *user_data, N_Vector tmp1,
                           N_Vector tmp2, N_Vector tmp3);

  
/*=================================================================
   IDASDLS Exported functions
  =================================================================*/

/*-----------------------------------------------------------------
  Required inputs for the IDASDLS linear solver interface:

  IDADlsSetLinearSolver specifies the direct SUNLinearSolver 
  object that IDAS should use.

  The return value is one of:
    IDADLS_SUCCESS   if successful
    IDADLS_MEM_NULL  if the IDA memory was NULL
    IDADLS_ILL_INPUT if the arguments are incompatible
---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDADlsSetLinearSolver(void *ida_mem, 
                                          SUNLinearSolver LS,
                                          SUNMatrix A);


/*---------------------------------------------------------------
  Optional inputs to the IDADLS linear solver interface:

  IDADlsSetJacFn specifies the dense Jacobian approximation
  routine to be used for a direct dense linear solver.
 
  By default, a difference quotient approximation is used for 
  dense and band; no default exists for sparse (so this must
  be user-supplied).
 
  The return value is one of:
     IDADLS_SUCCESS   if successful
     IDADLS_MEM_NULL  if the IDA memory was NULL
     IDADLS_LMEM_NULL if the linear solver memory was NULL
  -----------------------------------------------------------------*/
SUNDIALS_EXPORT int IDADlsSetJacFn(void *ida_mem, IDADlsJacFn jac);


/*---------------------------------------------------------------
 Optional outputs from the IDADLS linear solver

 IDADlsGetWorkSpace   returns the real and integer workspace used
                      by the direct linear solver.
 IDADlsGetNumJacEvals returns the number of calls made to the
                      Jacobian evaluation routine jac.
 IDADlsGetNumResEvals returns the number of calls to the user
                      F routine due to finite difference Jacobian
                      evaluation.
 IDADlsGetLastFlag    returns the last error flag set by any of
                      the IDADLS interface functions.

 The return value of IDADlsGet* is one of:
    IDADLS_SUCCESS   if successful
    IDADLS_MEM_NULL  if the IDA memory was NULL
    IDADLS_LMEM_NULL if the linear solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDADlsGetWorkSpace(void *ida_mem,
                                       long int *lenrwLS,
                                       long int *leniwLS);
SUNDIALS_EXPORT int IDADlsGetNumJacEvals(void *ida_mem,
                                         long int *njevals);
SUNDIALS_EXPORT int IDADlsGetNumResEvals(void *ida_mem,
                                         long int *nfevalsLS);
SUNDIALS_EXPORT int IDADlsGetLastFlag(void *ida_mem,
                                      long int *flag);

/*---------------------------------------------------------------
 The following function returns the name of the constant 
 associated with a IDADLS return flag
---------------------------------------------------------------*/
SUNDIALS_EXPORT char *IDADlsGetReturnFlagName(long int flag);

  
/*=================================================================
  PART II:  Backward Problems
  =================================================================*/

/*-----------------------------------------------------------------
  FUNCTION TYPES
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  Type: IDADlsJacFnB

  A Jacobian approximation function JacB for the adjoint (backward) 
  problem must have the prototype given below. 
  -----------------------------------------------------------------*/
typedef int (*IDADlsJacFnB)(realtype tt, realtype c_jB, N_Vector yy, 
                            N_Vector yp, N_Vector yyB, N_Vector ypB,
                            N_Vector rrB, SUNMatrix JacB,
                            void *user_dataB, N_Vector tmp1B,
                            N_Vector tmp2B, N_Vector tmp3B);


/*-----------------------------------------------------------------
  Type: IDADlsJacFnBS

  A Jacobian approximation function JacBS for the adjoint (backward) 
  problem, sensitivity-dependent case, must have the prototype given 
  below. 
  -----------------------------------------------------------------*/
typedef int (*IDADlsJacFnBS)(realtype tt, realtype c_jB, N_Vector yy, 
                             N_Vector yp, N_Vector *yS, N_Vector *ypS,
                             N_Vector yyB, N_Vector ypB, N_Vector rrB,
                             SUNMatrix JacB, void *user_dataB, 
                             N_Vector tmp1B, N_Vector tmp2B,
                             N_Vector tmp3B);


  
/*-----------------------------------------------------------------
  Exported functions
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  Required inputs for the IDASDLS linear solver interface:

  IDADlsSetLinearSolverB specifies the direct SUNLinearSolver 
  object that IDAS should use for backward integration. The 
  'which' argument is the int returned by IDACreateB.

  The return value is one of:
    IDADLS_SUCCESS   if successful
    IDADLS_MEM_NULL  if the IDA memory was NULL
    IDADLS_ILL_INPUT if the arguments are incompatible
---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDADlsSetLinearSolverB(void *ida_mem,
                                           int which,
                                           SUNLinearSolver LS,
                                           SUNMatrix A);

/*--------------------------------------------------------------------
  Functions: IDADlsSetJacFnB and IDADlsSetJacFnBS

  IDADlsSetJacFnB specifies the Jacobian function to be used by a
  IDASDLS linear solver for the backward integration phase, when
  the backward problem does not depend on forward sensitivities.
  IDADlsSetJacFnBS specifies the Jacobian function when the backward 
  problem does depend on sensitivities.
  The 'which' argument is the int returned by IDACreateB.
  --------------------------------------------------------------------*/
SUNDIALS_EXPORT int IDADlsSetJacFnB(void *ida_mem, int which,
                                    IDADlsJacFnB jacB);
SUNDIALS_EXPORT int IDADlsSetJacFnBS(void *ida_mem, int which,
                                     IDADlsJacFnBS jacBS);


#ifdef __cplusplus
}
#endif

#endif
