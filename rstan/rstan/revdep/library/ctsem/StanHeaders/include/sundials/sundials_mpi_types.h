/* -----------------------------------------------------------------
 * Programmer(s): Scott Cohen, Alan Hindmarsh, Radu Serban,
 *                Aaron Collier, and Slaven Peles @ LLNL
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
 * This header file contains definitions of MPI data types, which 
 * are used by MPI parallel vector implementations.
 * -----------------------------------------------------------------*/

#include <sundials/sundials_types.h>

/* define MPI data types */

#if defined(SUNDIALS_SINGLE_PRECISION)
  #define PVEC_REAL_MPI_TYPE MPI_FLOAT
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  #define PVEC_REAL_MPI_TYPE MPI_DOUBLE
#elif defined(SUNDIALS_EXTENDED_PRECISION)
  #define PVEC_REAL_MPI_TYPE MPI_LONG_DOUBLE
#endif

#if defined(SUNDIALS_INT64_T)
  #define PVEC_INTEGER_MPI_TYPE MPI_INT64_T
#elif defined(SUNDIALS_INT32_T)
  #define PVEC_INTEGER_MPI_TYPE MPI_INT32_T
#endif

