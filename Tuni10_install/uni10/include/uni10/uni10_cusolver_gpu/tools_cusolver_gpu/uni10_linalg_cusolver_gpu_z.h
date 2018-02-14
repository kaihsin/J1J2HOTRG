#ifndef __UNI10_LINALG_CUSOLVER_GPU_Z_H__
#define __UNI10_LINALG_CUSOLVER_GPU_Z_H__

#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_tools_cusolver_gpu.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_cusolver_gpu_dz.h"

namespace uni10{

  namespace uni10_linalg{

    // Blas
    //
    void vectorAdd_gpu(uni10_complex128 a, uni10_complex128* X, uni10_int incx, uni10_complex128* Y, uni10_int incy, uni10_uint64 N);   // Y = a*X + Y

    void vectorAdd_gpu(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N);// Y = Y + X

    void vectorSub_gpu(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N);// Y = Y - X

    void vectorMul_gpu(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N); // Y = Y * X, element-wise multiplication;

    void vectorScal_gpu(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);// X = a * X

    void vectorExp_gpu(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);

    uni10_complex128 vectorSum_gpu(uni10_complex128* X, uni10_uint64 N, uni10_int inc);

    uni10_double64 vectorNorm_gpu(uni10_complex128* X, uni10_uint64 N, uni10_int inc);

    void matrixDot_gpu(uni10_complex128* A, uni10_complex128* B, uni10_int M, uni10_int N, uni10_int K, uni10_complex128* C);

    void diagRowMul_gpu(uni10_complex128* mat, uni10_complex128* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul_gpu(uni10_complex128* mat, uni10_complex128* diag, uni10_uint64 M, uni10_uint64 N);

    void setTranspose_gpu(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N, uni10_complex128* AT);

    void setTranspose_gpu(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    void setDagger_gpu(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N, uni10_complex128* AT);

    void setDagger_gpu(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    void setConjugate_gpu(uni10_complex128 *A, uni10_uint64 N, uni10_complex128 *A_conj);

    void setConjugate_gpu(uni10_complex128 *A, uni10_uint64 N);

    void setIdentity_gpu(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    // Lapack
    //
    void matrixSVD_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* U, uni10_complex128* S, uni10_complex128* vT);

    void matrixSDD_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* U, uni10_complex128* S, uni10_complex128* vT);

    void matrixQR_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* R);

    void matrixRQ_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* R, uni10_complex128* Q);

    void matrixQL_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* L);

    void matrixLQ_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* L, uni10_complex128* Q);

    void matrixQDR_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* D, uni10_complex128* R);

    void matrixLDQ_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* L, uni10_complex128* D, uni10_complex128* Q);

    void matrixQDRCPIVOT_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* D, uni10_complex128* R);

    void matrixInv_gpu(uni10_complex128* A, uni10_int N);

    uni10_complex128 matrixDet_gpu(uni10_complex128* A, uni10_int N);

    //=================================================================================//

    void eigDecompose_gpu(uni10_complex128* Kij, uni10_int N, uni10_complex128* Eig, uni10_complex128 *EigVec);

    //
    // function overload for operators + - * += -= *=
    // 
    void matrix_diag_dense_add_gpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_sub_gpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_mul_gpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_add_gpu(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_dense_diag_sub_gpu(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_diag_dense_mul_gpu(std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n);
   
    void matrix_dense_diag_add_gpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_sub_gpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_mul_gpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
