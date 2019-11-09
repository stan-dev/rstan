/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer: Carol Woodward, Slaven Peles @ LLNL,
 *             Daniel R. Reynolds @ SMU.
 * -----------------------------------------------------------------
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This header file contains definitions and declarations for use by
 * sparse linear solvers for Ax = b. 
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_SPARSE_H
#define _SUNDIALS_SPARSE_H

#include <stdio.h>

#include <sundials/sundials_types.h>
#include <sundials/sundials_direct.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * ==================================================================
 * Type definitions
 * ==================================================================
 */

#define CSC_MAT 0
#define CSR_MAT 1

/*
 * -----------------------------------------------------------------
 * Type : SlsMat
 * -----------------------------------------------------------------
 * The type SlsMat is defined to be a pointer to a structure
 * with various sizes, a data field, and arrays for the row and 
 * column information for the sparse matrix entries.
 * The M and N fields indicates the number 
 * of rows and columns, respectively. The data field is a one 
 * dimensional array used for component storage. The NNZ field indicates
 * the number of nonzero entries in the matrix. The integer array, asub, 
 * holds the row index for each of the matrix entries.  The integer
 * array, xa, holds the index entry for the starting value of each column.
 * -----------------------------------------------------------------
 * The relevant fields in DlsMat are:
 *    M     - number of rows
 *    N     - number of columns
 *    NNZ   - the number of nonzero entries in the matrix
 *    NP    - number of index pointers
 *    data  - pointer to a contiguous block of realtype variables
 *    sparsetype - type of sparse matrix: compressed sparse column or row
 *    indexvals  - indices of each nonzero entry (columns or rows)
 *    indexptrs  - starting index of the first entry in data for each slice
 *    rowvals - pointer to row indices of each nonzero entry
 *    colptrs - pointer to starting indices in data array for each column
 *    colvals - pointer to column indices of each nonzero entry
 *    rowptrs - pointer to starting indices in data array for each row
 *
 * The nonzero entries of the matrix are stored in
 * compressed column format.  Row indices of entries in 
 * column j are stored in rowvals[colptrs[j]] through rowvals[colptrs[j+i]-1]
 * and corresponding numerical values of the matrix are stored 
 * in the same entries of data.
 * -----------------------------------------------------------------
 */

typedef struct _SlsMat {
  int M;
  int N;
  int NNZ;
  int NP;
  realtype *data;
  int sparsetype;
  int *indexvals;
  int *indexptrs;
  /* CSC indices */
  int **rowvals;
  int **colptrs;
  /* CSR indices */
  int **colvals;
  int **rowptrs;
} *SlsMat;

/*
 * ==================================================================
 * Exported function prototypes (functions working on SlsMat)
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function: SparseNewMat
 * -----------------------------------------------------------------
 * SparseNewMat allocates memory for a compressed column sparse
 * matrix with M rows, N columns, NNZ nonzeros and of sparsetype 
 * type (CSC or CSR matrix). SparseNewMat returns NULL if the 
 * request for matrix storage cannot be satisfied. See the above 
 * documentation for the type SlsMat for matrix storage details.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SlsMat SparseNewMat(int M, int N, int NNZ, int sparsetype);

/*
 * -----------------------------------------------------------------
 * Function: SparseFromDenseMat
 * -----------------------------------------------------------------
 * SlsConvertDense creates a new CSC matrix from an existing
 * dense/band matrix by copying all nonzero values into the sparse 
 * matrix structure.  SlsConvertDense returns NULL if the request 
 * for matrix storage cannot be satisfied. 
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SlsMat SparseFromDenseMat(const DlsMat A, int sparsetype);

/*
 * -----------------------------------------------------------------
 * Functions: SparseDestroyMat
 * -----------------------------------------------------------------
 * SparseDestroyMat frees the memory allocated by SparseNewMat
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SparseDestroyMat(SlsMat A);

/*
 * -----------------------------------------------------------------
 * Function : SparseSetMatToZero
 * -----------------------------------------------------------------
 * SetToZero sets all the elements of the sparse matrix A to 0.0.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SparseSetMatToZero(SlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: SparseCopyMat
 * -----------------------------------------------------------------
 * This function copies sparse matrix A into sparse matrix B.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SparseCopyMat(const SlsMat A, SlsMat B);

/*
 * -----------------------------------------------------------------
 * Functions: SparseScaleMat
 * -----------------------------------------------------------------
 * This function scales all data entries of a sparse matrix A 
 * by the realtype number in b.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SparseScaleMat(realtype b, SlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: SparseAddIdentityMat
 * -----------------------------------------------------------------
 * This function adds 1 to every diagonal entry of A.
 * Note that the resulting matrix may have more nonzero entries than 
 * the original.  This is accounted for, so that the return matrix 
 * may be larger than the one sent in.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SparseAddIdentityMat(SlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: SparseAddMat
 * -----------------------------------------------------------------
 * This function adds two sparse matrices: A = A+B.
 * Note that the resulting matrix may have more nonzero entries than
 * either of the original matrices.  This is accounted for, so that 
 * the return matrix may be larger than the ones sent in.  Upon 
 * successful completion, the return value is zero; otherwise 1 is 
 * returned.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SparseAddMat(SlsMat A, const SlsMat B);

/*
 * -----------------------------------------------------------------
 * Functions: SparseReallocMat
 * -----------------------------------------------------------------
 * This function reallocs internal arrays so that the resulting matrix 
 * holds colptrs[N] nonzeros.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SparseReallocMat(SlsMat A);

/*
 * -----------------------------------------------------------------
 * Functions: SparseMatvec
 * -----------------------------------------------------------------
 * This function computes the matrix-vector product, y=A*x, where A
 * is a sparse matrix of dimension MxN, x is a realtype array of 
 * length N, and y is a realtype array of length M. Upon successful
 * completion, the return value is zero; otherwise 1 is returned.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int SparseMatvec(const SlsMat A, const realtype *x, realtype *y);

/*
 * -----------------------------------------------------------------
 * Functions: SparsePrintMat
 * -----------------------------------------------------------------
 * This function prints the compressed column matrix information for 
 * matrix A to standard output.
 * It is intended as a debugging tool with small values of NNZ.
 * The elements are printed using the %g/%lg/%Lg option. 
 * A blank line is printed before and after the matrix.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void SparsePrintMat(const SlsMat A, FILE* outfile);




#ifdef __cplusplus
}
#endif

#endif
