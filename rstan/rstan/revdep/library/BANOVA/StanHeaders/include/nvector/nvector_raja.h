/* ----------------------------------------------------------------- 
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the header file for the RAJA implementation of the
 * NVECTOR module.
 *
 * Part I contains declarations specific to the RAJA
 * implementation of the supplied NVECTOR module.
 *
 * Part II contains the prototype for the constructor N_VNew_Raja
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
 *       N_VLinearSum_Raja(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_RAJA_H
#define _NVECTOR_RAJA_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

    
    
/*
 * -----------------------------------------------------------------
 * PART I: RAJA implementation of N_Vector
 * -----------------------------------------------------------------
 */

/* RAJA implementation of the N_Vector 'content' structure
   contains the length of the vector, a pointer to an array
   of 'realtype' components, and a flag indicating ownership of
   the data */

struct _N_VectorContent_Raja {};

typedef struct _N_VectorContent_Raja *N_VectorContent_Raja;



/*
 * -----------------------------------------------------------------
 * PART II: functions exported by nvector_raja
 * 
 * CONSTRUCTORS:
 *    N_VNew_Raja
 *    N_VNewEmpty_Raja
 *    N_VMake_Raja
 *    N_VCloneVectorArray_Raja
 *    N_VCloneVectorArrayEmpty_Raja
 * DESTRUCTORS:
 *    N_VDestroy_Raja
 *    N_VDestroyVectorArray_Raja
 * OTHER:
 *    N_VGetLength_Raja
 *    N_VGetHostArrayPointer_Raja
 *    N_VGetDeviceArrayPointer_Raja
 *    N_VPrint_Raja
 *    N_VPrintFile_Raja
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNew_Raja
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a RAJA vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_Raja(sunindextype vec_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_Raja
 * -----------------------------------------------------------------
 * This function creates a new RAJA N_Vector with an empty (NULL)
 * data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Raja(sunindextype vec_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_Raja
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a RAJA vector
 * with a user-supplied data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_Raja(N_VectorContent_Raja c);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArray_Raja
 * -----------------------------------------------------------------
 * This function creates an array of 'count' RAJA vectors by
 * cloning a given vector w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_Raja(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArrayEmpty_Raja
 * -----------------------------------------------------------------
 * This function creates an array of 'count' RAJA vectors each
 * with an empty (NULL) data array by cloning w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_Raja(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VDestroyVectorArray_Raja
 * -----------------------------------------------------------------
 * This function frees an array of RAJA vectors created with 
 * N_VCloneVectorArray_Raja or N_VCloneVectorArrayEmpty_Raja.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VDestroyVectorArray_Raja(N_Vector *vs, int count);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetLength_Raja
 * -----------------------------------------------------------------
 * This function returns the length of the vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunindextype N_VGetLength_Raja(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetHostArrayPointer_Raja
 * -----------------------------------------------------------------
 * This function returns pointer to the host raw data.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT realtype *N_VGetHostArrayPointer_Raja(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetDeviceArrayPointer_Raja
 * -----------------------------------------------------------------
 * This function returns pointer to the device raw data.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT realtype *N_VGetDeviceArrayPointer_Raja(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VCopyTotDevice_Raja
 * -----------------------------------------------------------------
 * This function copies host data to the device.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VCopyToDevice_Raja(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VCopyTotDevice_Raja
 * -----------------------------------------------------------------
 * This function copies vector data from the device to the host.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VCopyFromDevice_Raja(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_Raja
 * -----------------------------------------------------------------
 * This function prints the content of a RAJA vector to stdout.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrint_Raja(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrintFile_Raja
 * -----------------------------------------------------------------
 * This function prints the content of a RAJA vector to outfile.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrintFile_Raja(N_Vector v, FILE *outfile);

/*
 * -----------------------------------------------------------------
 * RAJA implementations of various useful vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_Raja(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Raja(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Raja(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Raja(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Raja(N_Vector v, sunindextype *lrw, sunindextype *liw);
SUNDIALS_EXPORT realtype *N_VGetArrayPointer_Raja(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_Raja(realtype *v_data, N_Vector v);
SUNDIALS_EXPORT void N_VLinearSum_Raja(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Raja(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Raja(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Raja(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Raja(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Raja(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Raja(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Raja(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Raja(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Raja(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Raja(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Raja(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Raja(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Raja(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Raja(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Raja(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Raja(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Raja(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Raja(N_Vector num, N_Vector denom);

#ifdef __cplusplus
}
#endif

#endif
