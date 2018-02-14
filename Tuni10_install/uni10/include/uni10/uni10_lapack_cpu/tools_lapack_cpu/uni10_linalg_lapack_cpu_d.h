#ifndef __UNI10_LINALG_LAPACK_CPU_D_H__
#define __UNI10_LINALG_LAPACK_CPU_D_H__

#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_tools_lapack_cpu.h"

#ifdef UNI_MKL
#include "mkl.h"
#else
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_lapack_wrapper_cpu.h"
#endif

namespace uni10{

  namespace uni10_linalg{

    // Blas 
    //
    void vectorAdd(uni10_double64 a, uni10_double64* X, uni10_int incx, uni10_double64* Y, uni10_int incy, uni10_uint64 N);   // Y = a*X + Y

    void vectorAdd(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y + X

    void vectorSub(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y - X

    void vectorMul(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y * X, element-wise multiplication;

    void vectorScal(uni10_double64 a, uni10_double64* X, uni10_uint64 N);  // X = a * X

    void vectorExp(uni10_double64 a, uni10_double64* X, uni10_uint64 N);

    uni10_double64 vectorSum(uni10_double64* X, uni10_uint64 N, uni10_int inc);

    uni10_double64 vectorNorm(uni10_double64* X, uni10_uint64 N, uni10_int inc);

    void matrixDot(uni10_double64* A, uni10_double64* B, uni10_int M, uni10_int N, uni10_int K, uni10_double64* C);

    void diagRowMul(uni10_double64* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul(uni10_double64* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    void setTranspose(uni10_double64* A, uni10_uint64 M, uni10_uint64 N, uni10_double64* AT);

    void setTranspose(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void setDagger(uni10_double64* A, uni10_uint64 M, uni10_uint64 N, uni10_double64 *AT);

    void setDagger(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void setIdentity(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void trimatrixEigH(uni10_double64* D, uni10_double64* E, uni10_int N, 
        uni10_double64* z=NULL, uni10_int LDZ=1);

    // Lapack
    //
    /*Generate a set of row vectors which form a othonormal basis
     *For the incoming matrix "elem", the number of row <= the number of column, M <= N
     */
    void matrixSVD(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* U, uni10_double64* S, uni10_double64* vT);

    void matrixSDD(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* U, uni10_double64* S, uni10_double64* vT);

    void matrixQR(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* R);

    void matrixRQ(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* R, uni10_double64* Q);

    void matrixQL(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* L);

    void matrixLQ(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* L, uni10_double64* Q);

    void matrixQDR(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* D, uni10_double64* R);

    void matrixLDQ(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* L, uni10_double64* D, uni10_double64* Q);

    void matrixQDRCPIVOT(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* D, uni10_double64* R);

    void matrixInv(uni10_double64* A, uni10_int N);

    uni10_double64 matrixDet(uni10_double64* A, uni10_int N);

    //=====================================================================================//
    //
    void reshapeElem(uni10_double64* elem, uni10_uint64* transOffset);

    void eigSyDecompose(uni10_double64* Kij, uni10_int N, uni10_double64* Eig, uni10_double64* EigVec);

    //
    // function overload for operator + - * += -= *= 
    //
    void matrix_diag_dense_add(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_diag_dense_sub(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_diag_dense_mul(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_dense_diag_add(double* a, const double* D, uni10_uint64 m, uni10_uint64 n);
                                                                                         
    void matrix_dense_diag_sub(double* a, const double* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_diag_dense_mul(double* D, const double* a, uni10_uint64 m, uni10_uint64 n);
                                                                                         
    void matrix_dense_diag_add(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_dense_diag_sub(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_dense_diag_mul(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
