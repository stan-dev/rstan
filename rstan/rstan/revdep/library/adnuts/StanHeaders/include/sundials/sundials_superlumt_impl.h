/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
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
 * Implementation header file for the SUNDIALS interface to the 
 * SuperLUMT linear solver.
 * -----------------------------------------------------------------
 */

#ifndef _SUNSLUMT_IMPL_H
#define _SUNSLUMT_IMPL_H

#ifndef _SLUMT_H
#define _SLUMT_H
/* #include "pdsp_defs.h" */
#include "slu_mt_ddefs.h"
#endif

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Definition of SLUMTData
 * -----------------------------------------------------------------
 */
 
typedef struct SLUMTDataRec {
 
  /* Structure for SuperLUMT-specific data */
 
  SuperMatrix *s_A, *s_AC, *s_L, *s_U, *s_B;
  Gstat_t *Gstat;
  int *perm_r, *perm_c;
  int num_threads;
  double diag_pivot_thresh; 
  superlumt_options_t *superlumt_options;

  int s_ordering;
 
} *SLUMTData;
 
#ifdef __cplusplus
} 
#endif 
 
#endif 
