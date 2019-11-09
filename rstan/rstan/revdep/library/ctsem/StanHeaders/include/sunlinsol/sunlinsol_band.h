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
 * This is the header file for the band implementation of the 
 * SUNLINSOL module.
 * 
 * Part I contains declarations specific to the band implementation
 * of the supplied SUNLINSOL module.
 * 
 * Part II contains the prototype for the constructor 
 * SUNBandLinearSolver as well as implementation-specific prototypes 
 * for various useful solver operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can 
 *     be found in the header file sundials_linearsolver.h.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_BAND_H
#define _SUNLINSOL_BAND_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_band.h>
#include <sunmatrix/sunmatrix_band.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: Band implementation of SUNLinearSolver
 *
 * The band implementation of the SUNLinearSolver 'content' 
 * structure contains:
 *     N -- size of the linear system
 *     pivots -- index array for partial pivoting in LU factorization
 *     last_flag -- last error return flag from internal setup/solve
 * -----------------------------------------------------------------
 */
  
struct _SUNLinearSolverContent_Band {
  sunindextype  N;
  sunindextype *pivots;
  long int last_flag;
};

typedef struct _SUNLinearSolverContent_Band *SUNLinearSolverContent_Band;

  
/*
 * -----------------------------------------------------------------
 * PART II: functions exported by sunlinsol_band
 * 
 * CONSTRUCTOR:
 *    SUNBandLinearSolver creates and allocates memory for a banded
 *    matrix solver
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNBandLinearSolver(N_Vector y,
                                                    SUNMatrix A);

/*
 * -----------------------------------------------------------------
 * band implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_Band(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_Band(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_Band(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_Band(SUNLinearSolver S, SUNMatrix A,
                                        N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_Band(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_Band(SUNLinearSolver S,
                                        long int *lenrwLS,
                                        long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_Band(SUNLinearSolver S);
  
#ifdef __cplusplus
}
#endif

#endif
