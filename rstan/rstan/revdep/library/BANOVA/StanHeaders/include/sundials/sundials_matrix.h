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
 * This is the header file for a generic matrix package.
 * It defines the SUNMatrix structure (_generic_SUNMatrix) which
 * contains the following fields:
 *   - an implementation-dependent 'content' field which contains
 *     the description and actual data of the matrix
 *   - an 'ops' filed which contains a structure listing operations
 *     acting on such matrices
 *
 * Part I of this file contains enumeration constants for all 
 * SUNDIALS-defined matrix types, as well as a generic type for 
 * user-supplied matrix types.
 *
 * Part II of this file contains type declarations for the
 * _generic_SUNMatrix and _generic_SUNMatrix_Ops structures, as well
 * as references to pointers to such structures (SUNMatrix).
 *
 * Part III of this file contains the prototypes for the matrix
 * functions which operate on SUNMatrix.
 *
 * At a minimum, a particular implementation of a SUNMatrix must
 * do the following:
 *  - specify the 'content' field of SUNMatrix,
 *  - implement the operations on those SUNMatrix,
 *  - provide a constructor routine for new SUNMatrix objects
 *
 * Additionally, a SUNMatrix implementation may provide the following:
 *  - macros to access the underlying SUNMatrix data
 *  - a routine to print the content of a SUNMatrix
 * -----------------------------------------------------------------
 */

#ifndef _SUNMATRIX_H
#define _SUNMATRIX_H
 
#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
 
  
/*
 * -----------------------------------------------------------------
 * I. Implemented SUNMatrix types
 * -----------------------------------------------------------------
 */
 
typedef enum {
  SUNMATRIX_DENSE, 
  SUNMATRIX_BAND, 
  SUNMATRIX_SPARSE, 
  SUNMATRIX_CUSTOM
} SUNMatrix_ID;

  
/*
 * -----------------------------------------------------------------
 * II. Generic definition of SUNMatrix 
 * -----------------------------------------------------------------
 */
 
/* Forward reference for pointer to SUNMatrix_Ops object */
typedef struct _generic_SUNMatrix_Ops *SUNMatrix_Ops;
 
/* Forward reference for pointer to SUNMatrix object */
typedef struct _generic_SUNMatrix *SUNMatrix;
 
/* Structure containing function pointers to matrix operations  */  
struct _generic_SUNMatrix_Ops {
  SUNMatrix_ID (*getid)(SUNMatrix);
  SUNMatrix    (*clone)(SUNMatrix);
  void         (*destroy)(SUNMatrix);
  int          (*zero)(SUNMatrix);
  int          (*copy)(SUNMatrix, SUNMatrix);
  int          (*scaleadd)(realtype, SUNMatrix, SUNMatrix);
  int          (*scaleaddi)(realtype, SUNMatrix);
  int          (*matvec)(SUNMatrix, N_Vector, N_Vector);
  int          (*space)(SUNMatrix, long int*, long int*);
};
 
/* A matrix is a structure with an implementation-dependent
   'content' field, and a pointer to a structure of matrix
   operations corresponding to that implementation.  */
struct _generic_SUNMatrix {
  void *content;
  struct _generic_SUNMatrix_Ops *ops;
};

  
/*
 * -----------------------------------------------------------------
 * III. Functions exported by SUNMatrix module
 *
 * SUNMatGetID
 *   Returns an identifier for the matrix type from the enumeration 
 *   SUNMatrix_ID.  This will be queried by a given linear solver to 
 *   assess compatibility.
 *
 * SUNMatClone
 *   Creates a new matrix of the same type as an existing matrix.
 *   It does not copy the matrix, but rather allocates storage for
 *   the new matrix.
 *
 * SUNMatDestroy
 *   Destroys a matrix created with SUNMatClone.
 *
 * SUNMatZero
 *   Sets all matrix entries to zero
 *
 * SUNMatScaleAdd
 *   Performs the operation A = c*A + B.  Returns an error if A 
 *   and B have different types and/or dimensions.
 *
 * SUNMatCopy
 *   Performs the operation B = A.  Should return an error if A and 
 *   B have different types and/or dimensions.
 *
 * SUNMatScaleAddI
 *   Performs the operation A = c*A + I.  Returns an error if A is 
 *   not a square matrix.
 *
 * SUNMatMatvec
 *   Performs the matrix-vector product y = A*x.  Returns an error if 
 *   A, x and/or y have incompatible types and/or dimensions, or if 
 *   x and y are the same vector.
 *
 * SUNMatSpace
 *   Returns the real and integer workspace requirement for the 
 *   SUNMatrix object.
 *
 * -----------------------------------------------------------------
 *
 * The following table lists the matrix functions used by
 * different modules in SUNDIALS. The symbols in the table
 * have the following meaning:
 *   D  -  called by dense linear solver modules
 *   B  -  called by band linear solver modules
 *   I  -  called by iterative linear solver modules
 *   S  -  called by sparse linear solver modules
 *   BP -  called by band preconditioner module
 *   DP -  called by band-block diagonal preconditioner module
 *
 *                                MODULES                  
 * MATRIX        -----------------------------------------------
 * FUNCTIONS     CVODE(S)      ARKode      IDA(S)       KINSOL    
 * ----------------------------------------------------------------
 *  GetID        D B S BP DP  D B S BP DP  D B S BP DP  D B S BP DP
 *  Clone        D B S BP DP  D B S BP DP  BP DP
 *  Destroy      D B S BP DP  D B S BP DP  BP DP
 *  Zero         D B S BP DP  D B S BP DP  D B S BP DP  D B BP DP
 *  Copy         D B S BP DP  D B S BP DP
 *  ScaleAddI    D B S BP DP  D B S BP DP
 *  ScaleAdd                  D B S BP DP
 *  Matvec*                   D B S
 *  Space        D B S BP DP  D B S BP DP  D B S BP DP  D B S BP DP
 * -----------------------------------------------------------------
 *  Note: MatrixMatvec is only called by ARKode when solving 
 *        problems having non-identity mass matrix
 * -----------------------------------------------------------------
 */
  
SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatClone(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAdd(realtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvec(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int SUNMatSpace(SUNMatrix A, long int *lenrw,
                                long int *leniw);
 
#ifdef __cplusplus
}
#endif
#endif
