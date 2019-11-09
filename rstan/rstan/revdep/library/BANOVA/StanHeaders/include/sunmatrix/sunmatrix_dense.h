/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * Based on code sundials_direct.h by: Radu Serban @ LLNL
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
 * This is the header file for the dense implementation of the 
 * SUNMATRIX module.
 * 
 * Part I contains declarations specific to the dense implementation
 * of the supplied SUNMATRIX module.
 * 
 * Part II defines accessor macros that allow the user to 
 * efficiently use this SUNMatrix type without making explicit
 * references to the underlying data structure.
 *
 * Part III contains the prototype for the constructor 
 * SUNDenseMatrix as well as implementation-specific prototypes 
 * for various useful matrix operations.
 *
 * Notes:
 *
 *   - The definition of the generic SUNMatrix structure can be found
 *     in the header file sundials_matrix.h.
 *
 *   - The definition of the type 'realtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype' and 'indextype'.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNMATRIX_DENSE_H
#define _SUNMATRIX_DENSE_H

#include <sundials/sundials_matrix.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: Dense implementation of SUNMatrix
 *
 * The dense implementation of the SUNMatrix 'content' structure
 * contains:
 *   M     - number of rows
 *   N     - number of columns
 *   data  - pointer to a contiguous block of realtype variables
 *   ldata - length of the data array = M*N
 *   cols  - array of pointers. cols[j] points to the first element 
 *           of the j-th column of the matrix in the array data.
 * The elements of a dense matrix are stored columnwise (i.e. columns 
 * are stored one on top of the other in memory); i.e. if A is a 
 * SUNMatrix_Dense object, then the (i,j)th element of A (with 
 * 0 <= i < M and 0 <= j < N) is given by (A->data)[j*M+i].
 *
 * -----------------------------------------------------------------
 */
  
struct _SUNMatrixContent_Dense {
  sunindextype M;
  sunindextype N;
  realtype *data;
  sunindextype ldata;
  realtype **cols;
};

typedef struct _SUNMatrixContent_Dense *SUNMatrixContent_Dense;

/*
 * -----------------------------------------------------------------
 * PART II: macros SM_CONTENT_D, SM_DATA_D, SM_ROWS_D, SM_COLUMNS_D, 
 *          SM_COLUMN_D, and SM_ELEMENT_D
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * SUNMatrix A;
 * SUNMatrixContent_Dense A_cont;
 * realtype *A_col_j, *A_data, **A_cols, A_ij;
 * sunindextype i, j, A_rows, A_columns, A_ldata;
 *
 * (1) SM_CONTENT_D
 *
 *     This macro gives access to the contents of the dense
 *     SUNMatrix
 *
 *     The assignment A_cont = SM_CONTENT_D(A) sets A_cont to be
 *     a pointer to the dense SUNMatrix content structure.
 *
 * (2) SM_DATA_D, SM_COLS_D, SM_LDATA_D, SM_ROWS_D, SM_COLUMNS_D
 *
 *     These macros give access to the individual parts of
 *     the content structure of a dense SUNMatrix.
 *
 *     The assignment A_data = SM_DATA_D(A) sets A_data to be
 *     a pointer to the first component of A. 
 *
 *     The assignment A_cols = SM_COLS_D(A) sets A_cols to be
 *     a pointer to the content's 'cols' entry.
 *
 *     The assignment A_ldata = SM_LDATA_D(A) sets A_ldata to be
 *     the length of the data array for A. 
 *
 *     The assignment A_rows = SM_ROWS_D(A) sets A_rows to be
 *     the number of rows in A.
 *
 *     The assignment A_columns = SM_COLUMNS_D(A) sets A_columns 
 *     to be the number of columns in A.
 *
 * (3) SM_COLUMN_D and SM_ELEMENT_D
 *
 *     These macros give access to the individual columns and 
 *     elements of a dense SUNMatrix, respectively.  In the 
 *     following, the entries of a SUNMatrix are indexed (i,j) 
 *     where i=0,...,M-1 and j=0,...,N-1.
 *     
 *     The assignment A_col_j = SM_COLUMN_D(A,j) sets A_col_j to 
 *     be a pointer to the jth column of the M-by-N dense
 *     matrix A, 0 <= j < N.  After the assignment, A_col_j may 
 *     be treated as an array indexed from 0 to M-1.
 *     The (i,j)-th element of A is thus referenced by col_j[i].
 *
 *     The assignment A_ij = SM_ELEMENT_D(A,i,j) sets A_ij to 
 *     the value of the (i,j)th element of the dense M-by-N matrix
 *     A, 0 <= i < M ; 0 <= j < N.  Similarly, the assignment 
 *     SM_ELEMENT_D(A,i,j) = A_ij sets the value of A_ij into the 
 *     (i,j) location of the matrix A.
 *
 * -----------------------------------------------------------------
 */

#define SM_CONTENT_D(A)     ( (SUNMatrixContent_Dense)(A->content) )

#define SM_ROWS_D(A)        ( SM_CONTENT_D(A)->M )

#define SM_COLUMNS_D(A)     ( SM_CONTENT_D(A)->N )

#define SM_LDATA_D(A)       ( SM_CONTENT_D(A)->ldata )

#define SM_DATA_D(A)        ( SM_CONTENT_D(A)->data )

#define SM_COLS_D(A)        ( SM_CONTENT_D(A)->cols )

#define SM_COLUMN_D(A,j)    ( (SM_CONTENT_D(A)->cols)[j] )

#define SM_ELEMENT_D(A,i,j) ( (SM_CONTENT_D(A)->cols)[j][i] )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by sunmatrix_dense
 * 
 * CONSTRUCTORS:
 *    SUNDenseMatrix
 * OTHER:
 *    SUNDenseMatrix_Print
 *    SUNDenseMatrix_Rows
 *    SUNDenseMatrix_Columns
 *    SUNDenseMatrix_LData
 *    SUNDenseMatrix_Data 
 *    SUNDenseMatrix_Cols
 *    SUNDenseMatrix_Column
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function: SUNDenseMatrix
 * -----------------------------------------------------------------
 * Creates and allocates memory for an M-by-N dense SUNMatrix.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNMatrix SUNDenseMatrix(sunindextype M, sunindextype N);

/*
 * -----------------------------------------------------------------
 * Functions: SUNDenseMatrix_Print
 * -----------------------------------------------------------------
 * This function prints the content of a M-by-N dense matrix A to
 * file pointer as it would normally appear on paper.
 * It is intended as debugging tools with small values of M and N.
 * The elements are printed using the %g/%lg/%Lg option. 
 * A blank line is printed before and after the matrix.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void SUNDenseMatrix_Print(SUNMatrix A, FILE* outfile);


/*
 * -----------------------------------------------------------------
 * Accessor Functions: 
 *
 * SUNDenseMatrix_Rows 
 *    Returns the number of rows in the dense matrix
 *
 * SUNDenseMatrix_Columns
 *    Returns the number of columns in the dense matrix
 *
 * SUNDenseMatrix_LData
 *    Returns the total allocated data length for the dense matrix
 *
 * SUNDenseMatrix_Data
 *    Returns a pointer to the data array for the dense matrix
 *
 * SUNDenseMatrix_Cols
 *    Returns a pointer to the cols array for the dense matrix
 *
 * SUNDenseMatrix_Column
 *    Returns a pointer to the jth column of the dense matrix
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunindextype SUNDenseMatrix_Rows(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNDenseMatrix_Columns(SUNMatrix A);
SUNDIALS_EXPORT sunindextype SUNDenseMatrix_LData(SUNMatrix A);
SUNDIALS_EXPORT realtype* SUNDenseMatrix_Data(SUNMatrix A);
SUNDIALS_EXPORT realtype** SUNDenseMatrix_Cols(SUNMatrix A);
SUNDIALS_EXPORT realtype* SUNDenseMatrix_Column(SUNMatrix A, sunindextype j);

/*
 * -----------------------------------------------------------------
 * dense implementations of various useful matrix operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNMatrix_ID SUNMatGetID_Dense(SUNMatrix A);
SUNDIALS_EXPORT SUNMatrix SUNMatClone_Dense(SUNMatrix A);
SUNDIALS_EXPORT void SUNMatDestroy_Dense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatZero_Dense(SUNMatrix A);
SUNDIALS_EXPORT int SUNMatCopy_Dense(SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAdd_Dense(realtype c, SUNMatrix A, SUNMatrix B);
SUNDIALS_EXPORT int SUNMatScaleAddI_Dense(realtype c, SUNMatrix A);
SUNDIALS_EXPORT int SUNMatMatvec_Dense(SUNMatrix A, N_Vector x, N_Vector y);
SUNDIALS_EXPORT int SUNMatSpace_Dense(SUNMatrix A, long int *lenrw,
                                      long int *leniw);

  
#ifdef __cplusplus
}
#endif

#endif
