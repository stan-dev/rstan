/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This header file is for routines to get SUNDIALS version info
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_VERSION_H
#define _SUNDIALS_VERSION_H

#include <sundials/sundials_config.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Function: SUNDIALSGetVersion
 * -----------------------------------------------------------------
 * SUNDIALSGetVersion fills a string with SUNDIALS version
 * information using the format major.minor.patch-label
 * 
 * Inputs:
 *   version = string to fill with version information
 *   len     = length of the input string version
 *
 * SUNDIALSGetVersion returns 0 if successful and -1 if the input
 * string is too short to store the SUNDIALS version.
 * -----------------------------------------------------------------*/
SUNDIALS_EXPORT int SUNDIALSGetVersion(char *version, int len);


/* -----------------------------------------------------------------
 * Function: SUNDIALSGetVersionNumber
 * -----------------------------------------------------------------
 * SUNDIALSGetVersionNumber fills separate integers with the 
 * SUNDIALS major, minor, and patch release version numbers and a
 * fills a string with the release label.
 * 
 * Inputs:
 *   major = integer for major version number
 *   minor = integer for minor version number
 *   patch = integer for patch version number
 *   label = string for version label
 *   len   = length of the input string label
 *
 * SUNDIALSGetVersionNumber returns 0 if successful and -1 if the
 * input string is too short to store the SUNDIALS version.
 * -----------------------------------------------------------------*/
SUNDIALS_EXPORT int SUNDIALSGetVersionNumber(int *major, int *minor, int *patch, 
                                             char *label, int len);

#ifdef __cplusplus
}
#endif

#endif
