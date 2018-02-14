#ifndef __UNI10_LINALG_LAPACK_CPU_Z_H__
#define __UNI10_LINALG_LAPACK_CPU_Z_H__

#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_tools_cusolver_gpu.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_lapack_cpu_dz.h"

#ifdef UNI_MKL
#include "mkl.h"
#else
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_lapack_wrapper_cpu.h"
#endif

namespace uni10{

  namespace uni10_linalg{

    // Blas
    //
    void vectorAdd_cpu(uni10_complex128 a, uni10_complex128* X, uni10_int incx, uni10_complex128* Y, uni10_int incy, uni10_uint64 N);   // Y = a*X + Y

    void vectorAdd_cpu(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N);// Y = Y + X

    void vectorSub_cpu(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N);// Y = Y - X

    void vectorMul_cpu(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N); // Y = Y * X, element-wise multiplication;

    void vectorScal_cpu(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);// X = a * X

    void vectorExp_cpu(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);

    uni10_complex128 vectorSum_cpu(uni10_complex128* X, uni10_uint64 N, uni10_int inc);

    uni10_double64 vectorNorm_cpu(uni10_complex128* X, uni10_uint64 N, uni10_int inc);

    void matrixDot_cpu(uni10_complex128* A, uni10_complex128* B, uni10_int M, uni10_int N, uni10_int K, uni10_complex128* C);

    void diagRowMul_cpu(uni10_complex128* mat, uni10_complex128* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul_cpu(uni10_complex128* mat, uni10_complex128* diag, uni10_uint64 M, uni10_uint64 N);

    void setTranspose_cpu(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N, uni10_complex128* AT);

    void setTranspose_cpu(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    void setDagger_cpu(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N, uni10_complex128* AT);

    void setDagger_cpu(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    void setConjugate_cpu(uni10_complex128 *A, uni10_uint64 N, uni10_complex128 *A_conj);

    void setConjugate_cpu(uni10_complex128 *A, uni10_uint64 N);

    void setIdentity_cpu(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    // Lapack
    //
    void matrixSVD_cpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* U, uni10_complex128* S, uni10_complex128* vT);

    void matrixSDD_cpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* U, uni10_complex128* S, uni10_complex128* vT);

    void matrixQR_cpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* R);

    void matrixRQ_cpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* R, uni10_complex128* Q);

    void matrixQL_cpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* L);

    void matrixLQ_cpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* L, uni10_complex128* Q);

    void matrixQDR_cpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* D, uni10_complex128* R);

    void matrixLDQ_cpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* L, uni10_complex128* D, uni10_complex128* Q);

    void matrixQDRCPIVOT_cpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* D, uni10_complex128* R);

    void matrixInv_cpu(uni10_complex128* A, uni10_int N);

    uni10_complex128 matrixDet_cpu(uni10_complex128* A, uni10_int N);

    //=================================================================================//

    void eigDecompose_cpu(uni10_complex128* Kij, uni10_int N, uni10_complex128* Eig, uni10_complex128 *EigVec);

    //
    // function overload for operators + - * += -= *=
    // 
    void matrix_diag_dense_add_cpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_sub_cpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_mul_cpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_add_cpu(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_dense_diag_sub_cpu(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_diag_dense_mul_cpu(std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n);
   
    void matrix_dense_diag_add_cpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_sub_cpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_mul_cpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);


  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
