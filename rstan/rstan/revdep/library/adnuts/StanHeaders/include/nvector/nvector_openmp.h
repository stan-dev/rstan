/* ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner and Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This NVECTOR module is based on the NVECTOR 
 *                   Serial module by Scott D. Cohen, Alan C. 
 *                   Hindmarsh, Radu Serban, and Aaron Collier 
 *                   @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the OpenMP implementation of the
 * NVECTOR module.
 *
 * Part I contains declarations specific to the OpenMP
 * implementation of the supplied NVECTOR module.
 *
 * Part II defines accessor macros that allow the user to
 * efficiently use the type N_Vector without making explicit
 * references to the underlying data structure.
 *
 * Part III contains the prototype for the constructor N_VNew_OpenMP
 * as well as implementation-specific prototypes for various useful
 * vector operations.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be found
 *     in the header file sundials_nvector.h.
 *
 *   - The definition of the type 'realtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype'.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *       N_VLinearSum_OpenMP(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_OPENMP_H
#define _NVECTOR_OPENMP_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: OPENMP implementation of N_Vector
 * -----------------------------------------------------------------
 */

/* OpenMP implementation of the N_Vector 'content' structure
   contains the length of the vector, a pointer to an array
   of 'realtype' components, and a flag indicating ownership of
   the data */

struct _N_VectorContent_OpenMP {
  sunindextype length;
  booleantype own_data;
  realtype *data;
  int num_threads;
};

typedef struct _N_VectorContent_OpenMP *N_VectorContent_OpenMP;

/*
 * -----------------------------------------------------------------
 * PART II: macros NV_CONTENT_OMP, NV_DATA_OMP, NV_OWN_DATA_OMP,
 *          NV_LENGTH_OMP, and NV_Ith_OMP
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * N_Vector v;
 * sunindextype i;
 *
 * (1) NV_CONTENT_OMP
 *
 *     This routines gives access to the contents of the OpenMP
 *     vector N_Vector.
 *
 *     The assignment v_cont = NV_CONTENT_OMP(v) sets v_cont to be
 *     a pointer to the OpenMP N_Vector content structure.
 *
 * (2) NV_DATA_OMP NV_OWN_DATA_OMP and NV_LENGTH_OMP
 *
 *     These routines give access to the individual parts of
 *     the content structure of a OpenMP N_Vector.
 *
 *     The assignment v_data = NV_DATA_OMP(v) sets v_data to be
 *     a pointer to the first component of v. The assignment
 *     NV_DATA_OMP(v) = data_V sets the component array of v to
 *     be data_v by storing the pointer data_v.
 *
 *     The assignment v_len = NV_LENGTH_OMP(v) sets v_len to be
 *     the length of v. The call NV_LENGTH_OMP(v) = len_v sets
 *     the length of v to be len_v.
 *
 * (3) NV_Ith_OMP
 *
 *     In the following description, the components of an
 *     N_Vector are numbered 0..n-1, where n is the length of v.
 *
 *     The assignment r = NV_Ith_OMP(v,i) sets r to be the value of
 *     the ith component of v. The assignment NV_Ith_OMP(v,i) = r
 *     sets the value of the ith component of v to be r.
 *
 * Note: When looping over the components of an N_Vector v, it is
 * more efficient to first obtain the component array via
 * v_data = NV_DATA_OMP(v) and then access v_data[i] within the
 * loop than it is to use NV_Ith_OMP(v,i) within the loop.
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_OMP(v)  ( (N_VectorContent_OpenMP)(v->content) )

#define NV_LENGTH_OMP(v)   ( NV_CONTENT_OMP(v)->length )

#define NV_NUM_THREADS_OMP(v)   ( NV_CONTENT_OMP(v)->num_threads )

#define NV_OWN_DATA_OMP(v) ( NV_CONTENT_OMP(v)->own_data )

#define NV_DATA_OMP(v)     ( NV_CONTENT_OMP(v)->data )

#define NV_Ith_OMP(v,i)    ( NV_DATA_OMP(v)[i] )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by nvector_OpenMP
 * 
 * CONSTRUCTORS:
 *    N_VNew_OpenMP
 *    N_VNewEmpty_OpenMP
 *    N_VMake_OpenMP
 *    N_VCloneVectorArray_OpenMP
 *    N_VCloneVectorArrayEmpty_OpenMP
 * DESTRUCTORS:
 *    N_VDestroy_OpenMP
 *    N_VDestroyVectorArray_OpenMP
 * OTHER:
 *    N_VGetLength_OpenMP
 *    N_VPrint_OpenMP
 *    N_VPrintFile_OpenMP
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNew_OpenMP
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a OpenMP vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_OpenMP(sunindextype vec_length, int num_threads);

/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_OpenMP
 * -----------------------------------------------------------------
 * This function creates a new OpenMP N_Vector with an empty (NULL)
 * data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_OpenMP(sunindextype vec_length, int num_threads);

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_OpenMP
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a OpenMP vector
 * with a user-supplied data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_OpenMP(sunindextype vec_length, realtype *v_data, int num_threads);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArray_OpenMP
 * -----------------------------------------------------------------
 * This function creates an array of 'count' OPENMP vectors by
 * cloning a given vector w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_OpenMP(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArrayEmpty_OpenMP
 * -----------------------------------------------------------------
 * This function creates an array of 'count' OPENMP vectors each
 * with an empty (NULL) data array by cloning w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_OpenMP(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VDestroyVectorArray_OpenMP
 * -----------------------------------------------------------------
 * This function frees an array of OPENMP vectors created with 
 * N_VCloneVectorArray_OpenMP or N_VCloneVectorArrayEmpty_OpenMP.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VDestroyVectorArray_OpenMP(N_Vector *vs, int count);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetLength_OpenMP
 * -----------------------------------------------------------------
 * This function returns number of vector elements.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunindextype N_VGetLength_OpenMP(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_OpenMP
 * -----------------------------------------------------------------
 * This function prints the content of a OpenMP vector to stdout.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrint_OpenMP(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrintFile_OpenMP
 * -----------------------------------------------------------------
 * This function prints the content of a OpenMP vector to outfile.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrintFile_OpenMP(N_Vector v, FILE *outfile);

/*
 * -----------------------------------------------------------------
 * OpenMP implementations of various useful vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_OpenMP(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_OpenMP(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_OpenMP(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_OpenMP(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_OpenMP(N_Vector v, sunindextype *lrw, sunindextype *liw);
SUNDIALS_EXPORT realtype *N_VGetArrayPointer_OpenMP(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_OpenMP(realtype *v_data, N_Vector v);
SUNDIALS_EXPORT void N_VLinearSum_OpenMP(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_OpenMP(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_OpenMP(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_OpenMP(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_OpenMP(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_OpenMP(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_OpenMP(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_OpenMP(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_OpenMP(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_OpenMP(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_OpenMP(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_OpenMP(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_OpenMP(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_OpenMP(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_OpenMP(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_OpenMP(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_OpenMP(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_OpenMP(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_OpenMP(N_Vector num, N_Vector denom);

#ifdef __cplusplus
}
#endif

#endif
