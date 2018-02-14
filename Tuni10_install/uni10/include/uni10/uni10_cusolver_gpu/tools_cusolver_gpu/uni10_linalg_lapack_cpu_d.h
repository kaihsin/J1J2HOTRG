#ifndef __UNI10_LINALG_LAPACK_CPU_D_H__
#define __UNI10_LINALG_LAPACK_CPU_D_H__

#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_tools_cusolver_gpu.h"

#ifdef UNI_MKL
#include "mkl.h"
#else
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_lapack_wrapper_cpu.h"
#endif

namespace uni10{

  namespace uni10_linalg{

    // Blas 
    //
    void vectorAdd_cpu(uni10_double64 a, uni10_double64* X, uni10_int incx, uni10_double64* Y, uni10_int incy, uni10_uint64 N);   // Y = a*X + Y

    void vectorAdd_cpu(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y + X

    void vectorSub_cpu(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y - X

    void vectorMul_cpu(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y * X, element-wise multiplication;

    void vectorScal_cpu(uni10_double64 a, uni10_double64* X, uni10_uint64 N);  // X = a * X

    void vectorExp_cpu(uni10_double64 a, uni10_double64* X, uni10_uint64 N);

    uni10_double64 vectorSum_cpu(uni10_double64* X, uni10_uint64 N, uni10_int inc);

    uni10_double64 vectorNorm_cpu(uni10_double64* X, uni10_uint64 N, uni10_int inc);

    void matrixDot_cpu(uni10_double64* A, uni10_double64* B, uni10_int M, uni10_int N, uni10_int K, uni10_double64* C);

    void diagRowMul_cpu(uni10_double64* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul_cpu(uni10_double64* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    void setTranspose_cpu(uni10_double64* A, uni10_uint64 M, uni10_uint64 N, uni10_double64* AT);

    void setTranspose_cpu(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void setDagger_cpu(uni10_double64* A, uni10_uint64 M, uni10_uint64 N, uni10_double64 *AT);

    void setDagger_cpu(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void setIdentity_cpu(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void trimatrixEigH_cpu(uni10_double64* D, uni10_double64* E, uni10_int N, 
        uni10_double64* z=NULL, uni10_int LDZ=1);

    // Lapack
    //
    /*Generate a set of row vectors which form a othonormal basis
     *For the incoming matrix "elem", the number of row <= the number of column, M <= N
     */
    void matrixSVD_cpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* U, uni10_double64* S, uni10_double64* vT);

    void matrixSDD_cpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* U, uni10_double64* S, uni10_double64* vT);

    void matrixQR_cpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* R);

    void matrixRQ_cpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* R, uni10_double64* Q);

    void matrixQL_cpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* L);

    void matrixLQ_cpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* L, uni10_double64* Q);

    void matrixQDR_cpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* D, uni10_double64* R);

    void matrixLDQ_cpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* L, uni10_double64* D, uni10_double64* Q);

    void matrixQDRCPIVOT_cpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* D, uni10_double64* R);

    void matrixInv_cpu(uni10_double64* A, uni10_int N);

    uni10_double64 matrixDet_cpu(uni10_double64* A, uni10_int N);

    //=====================================================================================//
    //
    void reshapeElem_cpu(uni10_double64* elem, uni10_uint64* transOffset);

    void eigSyDecompose_cpu(uni10_double64* Kij, uni10_int N, uni10_double64* Eig, uni10_double64* EigVec);

    //
    // function overload for operator + - * += -= *= 
    //
    void matrix_diag_dense_add_cpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_diag_dense_sub_cpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_diag_dense_mul_cpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_dense_diag_add_cpu(double* a, const double* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_dense_diag_sub_cpu(double* a, const double* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_diag_dense_mul_cpu(double* D, const double* a, uni10_uint64 m, uni10_uint64 n);

    void matrix_dense_diag_add_cpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_dense_diag_sub_cpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_dense_diag_mul_cpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
