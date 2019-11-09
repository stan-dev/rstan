/* -----------------------------------------------------------------
 * Programmer(s): Jean M. Sexton @ SMU
 *                Slaven Peles @ LLNL
 * ----------------------------------------------------------------- 
 * Based on work by: Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                   and Aaron Collier @ LLNL
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
 * This is the main header file for the MPI-enabled implementation
 * of the NVECTOR module.
 *
 * Part I contains declarations specific to the parallel
 * implementation of the supplied NVECTOR module.
 *
 * Part II contains the prototype for the constructor
 * N_VNew_ParHyp as well as implementation-specific prototypes
 * for various useful vector operations.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be
 *     found in the header file sundials_nvector.h.
 *
 *   - The definition of the type realtype can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type booleantype.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *        N_VLinearSum_ParHyp(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_PARHYP_H
#define _NVECTOR_PARHYP_H

#include <mpi.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_mpi_types.h>

/* hypre header files */
#include <_hypre_parcsr_mv.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*
 * -----------------------------------------------------------------
 * PART I: PARALLEL implementation of N_Vector               
 * -----------------------------------------------------------------
 */

/* 
 * Parallel implementation of the N_Vector 'content' structure
 * contains the global and local lengths of the vector, a pointer
 * to an array of 'realtype components', the MPI communicator,
 * and a flag indicating ownership of the data. 
 */
struct _N_VectorContent_ParHyp {
  sunindextype local_length;      /* local vector length         */
  sunindextype global_length;     /* global vector length        */
  booleantype own_parvector;  /* ownership of HYPRE vector   */
  MPI_Comm comm;              /* pointer to MPI communicator */

  hypre_ParVector *x;   /* The actual hypre_ParVector object */
};

typedef struct _N_VectorContent_ParHyp *N_VectorContent_ParHyp;


/*
 * -----------------------------------------------------------------
 * PART II: functions exported by nvector_ParHyp
 * 
 * CONSTRUCTORS:
 *    N_VNewEmpty_ParHyp
 *    N_VMake_ParHyp
 *    N_VCloneVectorArray_ParHyp
 *    N_VCloneVectorArrayEmpty_ParHyp
 * DESTRUCTORS:
 *    N_VDestroy_ParHyp
 *    N_VDestroyVectorArray_ParHyp
 * ACCESSOR FUNCTIONS:
 *    N_VGetVector_ParHyp
 * OTHER:
 *    N_VPrint_ParHyp
 *    N_VPrintFile_ParHyp
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_ParHyp
 * -----------------------------------------------------------------
 * This function creates a new hypre vector wrapper without the 
 * hypre vector itself.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_ParHyp(MPI_Comm comm, 
                                            sunindextype local_length,
                                            sunindextype global_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_ParHyp
 * -----------------------------------------------------------------
 * This function creates a hypre vector wrapper around user-supplied 
 * hypre vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_ParHyp(hypre_ParVector *x);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArray_ParHyp
 * -----------------------------------------------------------------
 * This function creates an array of 'count' N_Vectors by cloning a 
 * given vector w. Both, the wrapper and the underlying hypre vector
 * are cloned.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_ParHyp(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArrayEmpty_ParHyp
 * -----------------------------------------------------------------
 * This function creates an array of 'count' empty hypre vector
 * wrappers by cloning w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_ParHyp(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VDestroyVectorArray_ParHyp
 * -----------------------------------------------------------------
 * This function frees an array of N_Vector created with 
 * N_VCloneVectorArray_ParHyp or N_VCloneVectorArrayEmpty_ParHyp.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VDestroyVectorArray_ParHyp(N_Vector *vs, int count);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetVector_ParHyp
 * -----------------------------------------------------------------
 * Extracts underlying HYPRE vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT hypre_ParVector *N_VGetVector_ParHyp(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_ParHyp
 * -----------------------------------------------------------------
 * This function prints the local content of a parallel vector to
 * stdout.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrint_ParHyp(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrintFile_ParHyp
 * -----------------------------------------------------------------
 * This function prints the local content of a parallel vector to
 * outfile.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrintFile_ParHyp(N_Vector v, FILE *outfile);

/*
 * -----------------------------------------------------------------
 * parallel implementations of the vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_ParHyp(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_ParHyp(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_ParHyp(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_ParHyp(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_ParHyp(N_Vector v, sunindextype *lrw, sunindextype *liw);
SUNDIALS_EXPORT realtype *N_VGetArrayPointer_ParHyp(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_ParHyp(realtype *v_data, N_Vector v);
SUNDIALS_EXPORT void N_VLinearSum_ParHyp(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_ParHyp(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_ParHyp(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_ParHyp(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_ParHyp(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_ParHyp(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_ParHyp(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_ParHyp(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_ParHyp(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_ParHyp(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_ParHyp(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_ParHyp(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_ParHyp(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_ParHyp(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_ParHyp(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_ParHyp(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_ParHyp(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_ParHyp(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_ParHyp(N_Vector num, N_Vector denom);

#ifdef __cplusplus
}
#endif

#endif
