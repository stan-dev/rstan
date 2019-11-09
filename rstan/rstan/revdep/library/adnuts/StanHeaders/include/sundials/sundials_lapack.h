/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 *             Daniel Reynolds @ SMU
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
 * This is the header file for a generic package of direct matrix
 * operations for use with BLAS/LAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_LAPACK_H
#define _SUNDIALS_LAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * ==================================================================
 * Blas and Lapack functions
 * ==================================================================
 */

#if defined(SUNDIALS_F77_FUNC)

#define dcopy_f77       SUNDIALS_F77_FUNC(dcopy, DCOPY)
#define dscal_f77       SUNDIALS_F77_FUNC(dscal, DSCAL)
#define dgemv_f77       SUNDIALS_F77_FUNC(dgemv, DGEMV)
#define dtrsv_f77       SUNDIALS_F77_FUNC(dtrsv, DTRSV)
#define dsyrk_f77       SUNDIALS_F77_FUNC(dsyrk, DSKYR)

#define dgbtrf_f77      SUNDIALS_F77_FUNC(dgbtrf, DGBTRF)
#define dgbtrs_f77      SUNDIALS_F77_FUNC(dgbtrs, DGBTRS)
#define dgetrf_f77      SUNDIALS_F77_FUNC(dgetrf, DGETRF)
#define dgetrs_f77      SUNDIALS_F77_FUNC(dgetrs, DGETRS)
#define dgeqp3_f77      SUNDIALS_F77_FUNC(dgeqp3, DGEQP3)
#define dgeqrf_f77      SUNDIALS_F77_FUNC(dgeqrf, DGEQRF)
#define dormqr_f77      SUNDIALS_F77_FUNC(dormqr, DORMQR)
#define dpotrf_f77      SUNDIALS_F77_FUNC(dpotrf, DPOTRF)
#define dpotrs_f77      SUNDIALS_F77_FUNC(dpotrs, DPOTRS)

#define scopy_f77       SUNDIALS_F77_FUNC(scopy, SCOPY)
#define sscal_f77       SUNDIALS_F77_FUNC(sscal, SSCAL)
#define sgemv_f77       SUNDIALS_F77_FUNC(sgemv, SGEMV)
#define strsv_f77       SUNDIALS_F77_FUNC(strsv, STRSV)
#define ssyrk_f77       SUNDIALS_F77_FUNC(ssyrk, SSKYR)

#define sgbtrf_f77      SUNDIALS_F77_FUNC(sgbtrf, SGBTRF)
#define sgbtrs_f77      SUNDIALS_F77_FUNC(sgbtrs, SGBTRS)
#define sgetrf_f77      SUNDIALS_F77_FUNC(sgetrf, SGETRF)
#define sgetrs_f77      SUNDIALS_F77_FUNC(sgetrs, SGETRS)
#define sgeqp3_f77      SUNDIALS_F77_FUNC(sgeqp3, SGEQP3)
#define sgeqrf_f77      SUNDIALS_F77_FUNC(sgeqrf, SGEQRF)
#define sormqr_f77      SUNDIALS_F77_FUNC(sormqr, SORMQR)
#define spotrf_f77      SUNDIALS_F77_FUNC(spotrf, SPOTRF)
#define spotrs_f77      SUNDIALS_F77_FUNC(spotrs, SPOTRS)

#else

#define dcopy_f77       dcopy_
#define dscal_f77       dscal_
#define dgemv_f77       dgemv_
#define dtrsv_f77       dtrsv_
#define dsyrk_f77       dsyrk_

#define dgbtrf_f77      dgbtrf_
#define dgbtrs_f77      dgbtrs_
#define dgeqp3_f77      dgeqp3_
#define dgeqrf_f77      dgeqrf_
#define dgetrf_f77      dgetrf_
#define dgetrs_f77      dgetrs_
#define dormqr_f77      dormqr_
#define dpotrf_f77      dpotrf_
#define dpotrs_f77      dpotrs_

#define scopy_f77       scopy_
#define sscal_f77       sscal_
#define sgemv_f77       sgemv_
#define strsv_f77       strsv_
#define ssyrk_f77       ssyrk_

#define sgbtrf_f77      sgbtrf_
#define sgbtrs_f77      sgbtrs_
#define sgeqp3_f77      sgeqp3_
#define sgeqrf_f77      sgeqrf_
#define sgetrf_f77      sgetrf_
#define sgetrs_f77      sgetrs_
#define sormqr_f77      sormqr_
#define spotrf_f77      spotrf_
#define spotrs_f77      spotrs_
  
#endif

/* Level-1 BLAS */
  
extern void dcopy_f77(int *n, const double *x, const int *inc_x, double *y, const int *inc_y);
extern void dscal_f77(int *n, const double *alpha, double *x, const int *inc_x);

extern void scopy_f77(int *n, const float *x, const int *inc_x, float *y, const int *inc_y);
extern void sscal_f77(int *n, const float *alpha, float *x, const int *inc_x);
  
/* Level-2 BLAS */

extern void dgemv_f77(const char *trans, int *m, int *n, const double *alpha, const double *a, 
		      int *lda, const double *x, int *inc_x, const double *beta, double *y, int *inc_y, 
		      int len_trans);

extern void dtrsv_f77(const char *uplo, const char *trans, const char *diag, const int *n, 
		      const double *a, const int *lda, double *x, const int *inc_x, 
		      int len_uplo, int len_trans, int len_diag);

extern void sgemv_f77(const char *trans, int *m, int *n, const float *alpha, const float *a, 
		      int *lda, const float *x, int *inc_x, const float *beta, float *y, int *inc_y, 
		      int len_trans);

extern void strsv_f77(const char *uplo, const char *trans, const char *diag, const int *n, 
		      const float *a, const int *lda, float *x, const int *inc_x, 
		      int len_uplo, int len_trans, int len_diag);
  
/* Level-3 BLAS */

extern void dsyrk_f77(const char *uplo, const char *trans, const int *n, const int *k, 
		      const double *alpha, const double *a, const int *lda, const double *beta, 
		      const double *c, const int *ldc, int len_uplo, int len_trans);

extern void ssyrk_f77(const char *uplo, const char *trans, const int *n, const int *k, 
		      const float *alpha, const float *a, const int *lda, const float *beta, 
		      const float *c, const int *ldc, int len_uplo, int len_trans);
  
/* LAPACK */

extern void dgbtrf_f77(const int *m, const int *n, const int *kl, const int *ku, 
		       double *ab, int *ldab, int *ipiv, int *info);

extern void dgbtrs_f77(const char *trans, const int *n, const int *kl, const int *ku, const int *nrhs, 
		       double *ab, const int *ldab, int *ipiv, double *b, const int *ldb, 
		       int *info, int len_trans);


extern void dgeqp3_f77(const int *m, const int *n, double *a, const int *lda, int *jpvt, double *tau, 
		       double *work, const int *lwork, int *info);

extern void dgeqrf_f77(const int *m, const int *n, double *a, const int *lda, double *tau, double *work, 
		       const int *lwork, int *info);

extern void dgetrf_f77(const int *m, const int *n, double *a, int *lda, int *ipiv, int *info);

extern void dgetrs_f77(const char *trans, const int *n, const int *nrhs, double *a, const int *lda, 
		       int *ipiv, double *b, const int *ldb, int *info, int len_trans);


extern void dormqr_f77(const char *side, const char *trans, const int *m, const int *n, const int *k, 
		       double *a, const int *lda, double *tau, double *c, const int *ldc, 
		       double *work, const int *lwork, int *info, int len_side, int len_trans);

extern void dpotrf_f77(const char *uplo, const int *n, double *a, int *lda, int *info, int len_uplo);

extern void dpotrs_f77(const char *uplo, const int *n, const int *nrhs, double *a, const int *lda, 
		       double *b, const int *ldb, int * info, int len_uplo);


extern void sgbtrf_f77(const int *m, const int *n, const int *kl, const int *ku, 
		       float *ab, int *ldab, int *ipiv, int *info);

extern void sgbtrs_f77(const char *trans, const int *n, const int *kl, const int *ku, const int *nrhs, 
		       float *ab, const int *ldab, int *ipiv, float *b, const int *ldb, 
		       int *info, int len_trans);


extern void sgeqp3_f77(const int *m, const int *n, float *a, const int *lda, int *jpvt, float *tau, 
		       float *work, const int *lwork, int *info);

extern void sgeqrf_f77(const int *m, const int *n, float *a, const int *lda, float *tau, float *work, 
		       const int *lwork, int *info);

extern void sgetrf_f77(const int *m, const int *n, float *a, int *lda, int *ipiv, int *info);

extern void sgetrs_f77(const char *trans, const int *n, const int *nrhs, float *a, const int *lda, 
		       int *ipiv, float *b, const int *ldb, int *info, int len_trans);


extern void sormqr_f77(const char *side, const char *trans, const int *m, const int *n, const int *k, 
		       float *a, const int *lda, float *tau, float *c, const int *ldc, 
		       float *work, const int *lwork, int *info, int len_side, int len_trans);

extern void spotrf_f77(const char *uplo, const int *n, float *a, int *lda, int *info, int len_uplo);

extern void spotrs_f77(const char *uplo, const int *n, const int *nrhs, float *a, const int *lda, 
		       float *b, const int *ldb, int * info, int len_uplo);


#ifdef __cplusplus
}
#endif

#endif
