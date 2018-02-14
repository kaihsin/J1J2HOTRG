#ifndef __UNI10_LINALG_CUSOLVER_GPU_DZ_H__
#define __UNI10_LINALG_CUSOLVER_GPU_DZ_H__

#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_tools_cusolver_gpu.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_cusolver_gpu_z.h"

namespace uni10{

  namespace uni10_linalg{

    //Blas
    //
    void vectorAdd_gpu(uni10_double64 a, uni10_double64* X, uni10_int incx, uni10_complex128* Y, uni10_int incy, uni10_uint64 N);   // Y = a*X + Y

    void vectorAdd_gpu(uni10_complex128* Y, uni10_double64* X, uni10_uint64 N);

    void vectorSub_gpu(uni10_complex128* Y, uni10_double64* X, uni10_uint64 N);

    void vectorMul_gpu(uni10_complex128* Y, uni10_double64* X, uni10_uint64 N);
    
    void vectorScal_gpu(uni10_double64 a, uni10_complex128* X, uni10_uint64 N);

    void vectorExp_gpu(uni10_double64 a, uni10_complex128* X, uni10_uint64 N);
 
    void matrixDot_gpu(uni10_double64* A, uni10_complex128* B, uni10_int M, uni10_int N, uni10_int K, uni10_complex128* C);

    void matrixDot_gpu(uni10_complex128* A, uni10_double64* B, uni10_int M, uni10_int N, uni10_int K, uni10_complex128* C);

    void diagRowMul_gpu(uni10_complex128* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul_gpu(uni10_complex128* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    //Lapack
    //

    //====================================================================//

    void eigDecompose_gpu(uni10_double64* Kij, uni10_int N, uni10_complex128* Eig, uni10_complex128 *EigVec);

    void eigSyDecompose_gpu(uni10_complex128* Kij, uni10_int N, uni10_double64* Eig, uni10_complex128* EigVec);

    void matrixSVD_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* U, uni10_double64 *S, uni10_complex128* vT);

    void matrixSDD_gpu(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* U, uni10_double64 *S, uni10_complex128* vT);

    //
    // function overload for operators + - * += -= *=
    // 
    // r operator(+,-,*,+=,-=,*=) z
    void matrix_diag_dense_add_gpu(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_sub_gpu(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_mul_gpu(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);
   
    void matrix_dense_diag_add_gpu(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_sub_gpu(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_mul_gpu(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    // z operator(+,-,*,+=,-=,*=) r

    void matrix_diag_dense_add_gpu(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_sub_gpu(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_mul_gpu(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_add_gpu(std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_dense_diag_sub_gpu(std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_diag_dense_mul_gpu(std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n);
   
    void matrix_dense_diag_add_gpu(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_sub_gpu(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_mul_gpu(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
