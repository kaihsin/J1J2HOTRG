/****************************************************************************
*  @file uni10_lapack_wrapper.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2016
*    National Taiwan University
*    National Tsing-Hua University
*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief C wrapper functions for fortran BLAS and LAPACK libraries
*  @author Ying-Jer Kao, Yun-Hsuan Chou
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef __UNI10_SCALAPACK_WRAPPER_MPI_H__
#define __UNI10_SCALAPACK_WRAPPER_MPI_H__

#include <stdint.h>

#include <complex>

#if defined(UNI_MKL)
#include "mkl_blas.h"
#include "mkl_pblas.h"
#include "mkl_scalapack.h"
#else

#include "uni10/uni10_type.h"

extern "C"{

void   daxpy_(const uni10_int *n, const double *alpha, const double *x, 
    const uni10_int *incx, double *y, const uni10_int *incy);
void   zaxpy_(const uni10_int *n, const uni10_complex128 *alpha, const uni10_complex128 *x, 
    const uni10_int *incx, uni10_complex128 *y, const uni10_int *incy);

void   dscal_(const uni10_int *n, const double *a, double *x, const uni10_int *incx);
void   zscal_(const uni10_int *n, const uni10_complex128 *a, uni10_complex128 *x, const uni10_int *incx);
void   zdscal_(const uni10_int *n, const double *a, uni10_complex128 *x, const uni10_int *incx);

double pdlange_(const char* norm, const uni10_int* m, const uni10_int* n, const double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* work);
double pzlange_(const char* norm, const uni10_int* m, const uni10_int* n, const uni10_complex128* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* work);

void pdasum_( uni10_int *n, double *asum, double *x, uni10_int *ix, uni10_int *jx, uni10_int *descx, uni10_int *incx );

void pdscal_( uni10_int *n, double *a, double *x, uni10_int *ix, uni10_int *jx, uni10_int *descx, uni10_int *incx );
void pzscal( uni10_int *n, double *a, double *x, uni10_int *ix, uni10_int *jx, uni10_int *descx, uni10_int *incx );

void pdgeqrf_(const uni10_int* m, const uni10_int* n, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* tau, double* work, const uni10_int* lwork, uni10_int* info);
void pdgeqlf_(const uni10_int* m, const uni10_int* n, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* tau, double* work, const uni10_int* lwork, uni10_int* info);
void pdgerqf_(const uni10_int* m, const uni10_int* n, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* tau, double* work, const uni10_int* lwork, uni10_int* info);
void pdgelqf_(const uni10_int* m, const uni10_int* n, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* tau, double* work, const uni10_int* lwork, uni10_int* info);

void	pdorgqr_(const uni10_int* m, const uni10_int* n, const uni10_int* k, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, const double* tau, double* work, const uni10_int* lwork, uni10_int* info);
void	pdorgrq_(const uni10_int* m, const uni10_int* n, const uni10_int* k, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, const double* tau, double* work, const uni10_int* lwork, uni10_int* info);
void	pdorgql_(const uni10_int* m, const uni10_int* n, const uni10_int* k, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, const double* tau, double* work, const uni10_int* lwork, uni10_int* info);
void	pdorglq_(const uni10_int* m, const uni10_int* n, const uni10_int* k, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, const double* tau, double* work, const uni10_int* lwork, uni10_int* info);

void pdlacpy_(const char* uplo, const uni10_int* m, const uni10_int* n, const double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* b, const uni10_int* ib, const uni10_int* jb, const uni10_int* descb);

void pdgeadd_( const char *trans, const int *m, const int *n, const double *alpha, 
    const double *a, const int *ia, const int *ja, const int *desca, const double *beta, 
    double *c, const int *ic, const int *jc, const int *descc  );
void pzgeadd_( const char *trans, const int *m, const int *n, const double *alpha, 
    const double *a, const int *ia, const int *ja, const int *desca, const double *beta, 
    double *c, const int *ic, const int *jc, const int *descc  );

void pdgemm_( const char *transa, const char *transb, const uni10_int *m, const uni10_int *n, const uni10_int *k, 
    const double *alpha, const double *a, const uni10_int *ia, const uni10_int *ja, const uni10_int *desca, 
    const double *b, const uni10_int *ib, const uni10_int *jb, const uni10_int *descb, 
    const double *beta, double *c, const uni10_int *ic, const uni10_int *jc, const uni10_int *descc );
void pzgemm_( const char *transa, const char *transb, const uni10_int *m, const uni10_int *n, const uni10_int *k, 
    const uni10_complex128 *alpha, const uni10_complex128 *a, const uni10_int *ia, const uni10_int *ja, const uni10_int *desca, 
    const uni10_complex128 *b, const uni10_int *ib, const uni10_int *jb, const uni10_int *descb, 
    const uni10_complex128 *beta, uni10_complex128 *c, const uni10_int *ic, const uni10_int *jc, const uni10_int *descc );

void pdgesv_(const uni10_int* n, const uni10_int* nrhs, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, 
    uni10_int* ipiv, double* b, const uni10_int* ib, const uni10_int* jb, const uni10_int* descb, uni10_int* info);
void pzgesv_(const uni10_int* n, const uni10_int* nrhs, uni10_complex128* a, const uni10_int* ia, const uni10_int* ja, 
    const uni10_int* desca, uni10_int* ipiv, uni10_complex128* b, const uni10_int* ib, const uni10_int* jb, const uni10_int* descb, uni10_int* info);

void  pdgesvd_(const char* jobu, const char* jobvt, const uni10_int* m, const uni10_int* n, const double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* s, double* u, const uni10_int* iu, const uni10_int* ju, const uni10_int* descu, double* vt, const uni10_int* ivt, const uni10_int* jvt, const uni10_int* descvt, double* work, const uni10_int* lwork, uni10_int* info);
void  pzgesvd_(const char* jobu, const char* jobvt, const uni10_int* m, const uni10_int* n, const uni10_complex128* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* s, uni10_complex128* u, const uni10_int* iu, const uni10_int* ju, const uni10_int* descu, uni10_complex128* vt, const uni10_int* ivt, const uni10_int* jvt, const uni10_int* descvt, uni10_complex128* work, const uni10_int* lwork, double* rwork, uni10_int* info);

};


inline void daxpy(const uni10_int *n, const double *alpha, const double *x, const uni10_int *incx, double *y, const uni10_int *incy)
{
  daxpy_(n, alpha, x, incx, y, incy);
}

inline void zaxpy(const uni10_int *n, const uni10_complex128 *alpha, const uni10_complex128 *x, const uni10_int *incx, uni10_complex128 *y, const uni10_int *incy)
{
  zaxpy_(n, alpha, x, incx, y, incy);
}

inline void dscal(const uni10_int *n, const double *a, double *x, const uni10_int *incx)
{
  dscal_(n, a, x, incx);
}
inline void zscal(const uni10_int *n, const uni10_complex128 *a, uni10_complex128 *x, const uni10_int *incx)
{
  zscal_(n, a, x, incx);
}

inline void zdscal(const uni10_int *n, const double *a, uni10_complex128 *x, const uni10_int *incx)
{
  zdscal_(n, a, x, incx);
}

inline double pdlange(const char* norm, const uni10_int* m, const uni10_int* n, const double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* work)
{
  return pdlange_(norm, m, n, a, ia, ja, desca, work);
}

inline double pzlange(const char* norm, const uni10_int* m, const uni10_int* n, const uni10_complex128* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* work)
{
  return pzlange_(norm, m, n, a, ia, ja, desca, work);
}

inline void pdasum( uni10_int *n, double *asum, double *x, uni10_int *ix, uni10_int *jx, uni10_int *descx, uni10_int *incx )
{
  pdasum_(n, asum, x, ix, jx, descx, incx);
}

inline void pdgemm( const char *transa, const char *transb, const uni10_int *m, const uni10_int *n, const uni10_int *k, const double *alpha, const double *a, const uni10_int *ia, const uni10_int *ja, const uni10_int *desca, const double *b, const uni10_int *ib, const uni10_int *jb, const uni10_int *descb, const double *beta, double *c, const uni10_int *ic, const uni10_int *jc, const uni10_int *descc )
{
  pdgemm_(transa, transb, m, n, k , alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc);
}

void pdgesv(const uni10_int* n, const uni10_int* nrhs, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, 
    uni10_int* ipiv, double* b, const uni10_int* ib, const uni10_int* jb, const uni10_int* descb, uni10_int* info)
{
  pdgesv_(n, nrhs, a, ia, ja, desca, ipiv, b, ib, jb, descb, info);
}
void pzgesv(const uni10_int* n, const uni10_int* nrhs, uni10_complex128* a, const uni10_int* ia, const uni10_int* ja, 
    const uni10_int* desca, uni10_int* ipiv, uni10_complex128* b, const uni10_int* ib, const uni10_int* jb, const uni10_int* descb, uni10_int* info)
{
  pzgesv_(n, nrhs, a, ia, ja, desca, ipiv, b, ib, jb, descb, info);
}

inline void pzgemm( const char *transa, const char *transb, const uni10_int *m, const uni10_int *n, const uni10_int *k, const uni10_complex128 *alpha, const uni10_complex128 *a, const uni10_int *ia, const uni10_int *ja, const uni10_int *desca, const uni10_complex128 *b, const uni10_int *ib, const uni10_int *jb, const uni10_int *descb, const uni10_complex128 *beta, uni10_complex128 *c, const uni10_int *ic, const uni10_int *jc, const uni10_int *descc )
{
  pzgemm_(transa, transb, m, n, k , alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc);
}

inline void pdgesvd(const char* jobu, const char* jobvt, const uni10_int* m, const uni10_int* n, const double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* s, double* u, const uni10_int* iu, const uni10_int* ju, const uni10_int* descu, double* vt, const uni10_int* ivt, const uni10_int* jvt, const uni10_int* descvt, double* work, const uni10_int* lwork, uni10_int* info)
{
  pdgesvd_(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, work, lwork, info);
}

inline void pzgesvd(const char* jobu, const char* jobvt, const uni10_int* m, const uni10_int* n, const uni10_complex128* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* s, uni10_complex128* u, const uni10_int* iu, const uni10_int* ju, const uni10_int* descu, uni10_complex128* vt, const uni10_int* ivt, const uni10_int* jvt, const uni10_int* descvt, uni10_complex128* work, const uni10_int* lwork, double* rwork, uni10_int* info)
{

  pzgesvd_(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, work, lwork, rwork, info);

}

inline void pdgeadd( const char *trans, const int *m, const int *n, const double *alpha, 
    const double *a, const int *ia, const int *ja, const int *desca, const double *beta, 
    double *c, const int *ic, const int *jc, const int *descc  )
{
  pdgeadd_(trans, m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc);
}

inline void pzgeadd( const char *trans, const int *m, const int *n, const double *alpha, 
    const double *a, const int *ia, const int *ja, const int *desca, const double *beta, 
    double *c, const int *ic, const int *jc, const int *descc  )
{
  pzgeadd_(trans, m, n, alpha, a, ia, ja, desca, beta, c, ic, jc, descc);
}

inline void pdgeqrf(const uni10_int* m, const uni10_int* n, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* tau, double* work, const uni10_int* lwork, uni10_int* info)
{
  pdgeqrf_(m, n, a, ia, ja, desca, tau, work, lwork, info);
}

inline void pdgeqlf(const uni10_int* m, const uni10_int* n, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* tau, double* work, const uni10_int* lwork, uni10_int* info)
{
  pdgeqlf_(m, n, a, ia, ja, desca, tau, work, lwork, info);

}

inline void pdgerqf(const uni10_int* m, const uni10_int* n, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* tau, double* work, const uni10_int* lwork, uni10_int* info)
{
  pdgerqf_(m, n, a, ia, ja, desca, tau, work, lwork, info);
}

inline void pdgelqf(const uni10_int* m, const uni10_int* n, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* tau, double* work, const uni10_int* lwork, uni10_int* info)
{
  pdgelqf_(m, n, a, ia, ja, desca, tau, work, lwork, info);
}

inline void pdorgqr(const uni10_int* m, const uni10_int* n, const uni10_int* k, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, const double* tau, double* work, const uni10_int* lwork, uni10_int* info)
{
  pdorgqr_(m, n, k, a, ia, ja, desca, tau, work, lwork, info);

}

inline void pdorgrq(const uni10_int* m, const uni10_int* n, const uni10_int* k, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, const double* tau, double* work, const uni10_int* lwork, uni10_int* info)
{
  pdorgrq_(m, n, k, a, ia, ja, desca, tau, work, lwork, info);
}

inline void pdorgql(const uni10_int* m, const uni10_int* n, const uni10_int* k, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, const double* tau, double* work, const uni10_int* lwork, uni10_int* info)
{
  pdorgql_(m, n, k, a, ia, ja, desca, tau, work, lwork, info);
}

inline void pdorglq(const uni10_int* m, const uni10_int* n, const uni10_int* k, double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, const double* tau, double* work, const uni10_int* lwork, uni10_int* info)
{
  pdorglq_(m, n, k, a, ia, ja, desca, tau, work, lwork, info);
}

inline void pdlacpy(const char* uplo, const uni10_int* m, const uni10_int* n, const double* a, const uni10_int* ia, const uni10_int* ja, const uni10_int* desca, double* b, const uni10_int* ib, const uni10_int* jb, const uni10_int* descb)
{
  pdlacpy_(uplo, m, n, a, ia, ja, desca, b, ib, jb, descb);
}

#endif // End of wrapper.

#endif

