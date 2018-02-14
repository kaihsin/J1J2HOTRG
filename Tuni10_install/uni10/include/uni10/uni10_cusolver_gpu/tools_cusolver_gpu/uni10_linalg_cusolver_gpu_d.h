#ifndef __UNI10_LINALG_CUSOLVER_GPU_D_H__
#define __UNI10_LINALG_CUSOLVER_GPU_D_H__

#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_tools_cusolver_gpu.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/cuda_kernel_funcs/uni10_kernel_gpu.h"

namespace uni10{

  namespace uni10_linalg{

    // Blas 
    //
    void vectorAdd_gpu(uni10_double64 a, uni10_double64* X, uni10_int incx, uni10_double64* Y, uni10_int incy, uni10_uint64 N);   // Y = a*X + Y

    void vectorAdd_gpu(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y + X

    void vectorSub_gpu(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y - X

    void vectorMul_gpu(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y * X, element-wise multiplication;

    void vectorScal_gpu(uni10_double64 a, uni10_double64* X, uni10_uint64 N);  // X = a * X

    void vectorExp_gpu(uni10_double64 a, uni10_double64* X, uni10_uint64 N);

    uni10_double64 vectorSum_gpu(uni10_double64* X, uni10_uint64 N, uni10_int inc);

    uni10_double64 vectorNorm_gpu(uni10_double64* X, uni10_uint64 N, uni10_int inc);

    void matrixDot_gpu(uni10_double64* A, uni10_double64* B, uni10_int M, uni10_int N, uni10_int K, uni10_double64* C);

    void diagRowMul_gpu(uni10_double64* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul_gpu(uni10_double64* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    void setTranspose_gpu(uni10_double64* A, uni10_uint64 M, uni10_uint64 N, uni10_double64* AT);

    void setTranspose_gpu(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void setDagger_gpu(uni10_double64* A, uni10_uint64 M, uni10_uint64 N, uni10_double64 *AT);

    void setDagger_gpu(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void setIdentity_gpu(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void trimatrixEigH_gpu(uni10_double64* D, uni10_double64* E, uni10_int N, 
        uni10_double64* z=NULL, uni10_int LDZ=1);

    // Lapack
    //
    /*Generate a set of row vectors which form a othonormal basis
     *For the incoming matrix "elem", the number of row <= the number of column, M <= N
     */
    void matrixSVD_gpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* U, uni10_double64* S, uni10_double64* vT);

    void matrixSDD_gpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* U, uni10_double64* S, uni10_double64* vT);

    void matrixQR_gpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* R);

    void matrixRQ_gpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* R, uni10_double64* Q);

    void matrixQL_gpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* L);

    void matrixLQ_gpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* L, uni10_double64* Q);

    void matrixQDR_gpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* D, uni10_double64* R);

    void matrixLDQ_gpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* L, uni10_double64* D, uni10_double64* Q);

    void matrixQDRCPIVOT_gpu(uni10_double64* Mij_ori, uni10_int M, uni10_int N, uni10_double64* Q, uni10_double64* D, uni10_double64* R);

    void matrixInv_gpu(uni10_double64* A, uni10_int N);

    uni10_double64 matrixDet_gpu(uni10_double64* A, uni10_int N);

    //=====================================================================================//
    //
    void reshapeElem_gpu(uni10_double64* elem, uni10_uint64* transOffset);

    void eigSyDecompose_gpu(uni10_double64* Kij, uni10_int N, uni10_double64* Eig, uni10_double64* EigVec);

    //
    // function overload for operator + - * += -= *= 
    //
    void matrix_diag_dense_add_gpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_diag_dense_sub_gpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_diag_dense_mul_gpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_dense_diag_add_gpu(double* a, const double* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_dense_diag_sub_gpu(double* a, const double* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_diag_dense_mul_gpu(double* D, const double* a, uni10_uint64 m, uni10_uint64 n);

    void matrix_dense_diag_add_gpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_dense_diag_sub_gpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b);

    void matrix_dense_diag_mul_gpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
