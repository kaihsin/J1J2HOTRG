#ifndef __UNI10_LINALG_HYBRID_Z_H__
#define __UNI10_LINALG_HYBRID_Z_H__

#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_lapack_cpu_z.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_cusolver_gpu_z.h"

namespace uni10{

  namespace uni10_linalg{

    // Blas
    //
    void vectorAdd(uni10_complex128 a, uni10_complex128* X, uni10_int incx, uni10_complex128* Y, uni10_int incy, uni10_uint64 N);   // Y = a*X + Y

    void vectorAdd(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N);// Y = Y + X

    void vectorSub(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N);// Y = Y - X

    void vectorMul(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N); // Y = Y * X, element-wise multiplication;

    void vectorScal(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);// X = a * X

    void vectorExp(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);

    uni10_complex128 vectorSum(uni10_complex128* X, uni10_uint64 N, uni10_int inc);

    uni10_double64 vectorNorm(uni10_complex128* X, uni10_uint64 N, uni10_int inc);

    void matrixDot(uni10_complex128* A, uni10_complex128* B, uni10_int M, uni10_int N, uni10_int K, uni10_complex128* C);

    void diagRowMul(uni10_complex128* mat, uni10_complex128* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul(uni10_complex128* mat, uni10_complex128* diag, uni10_uint64 M, uni10_uint64 N);

    void setTranspose(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N, uni10_complex128* AT);

    void setTranspose(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    void setDagger(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N, uni10_complex128* AT);

    void setDagger(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    void setConjugate(uni10_complex128 *A, uni10_uint64 N, uni10_complex128 *A_conj);

    void setConjugate(uni10_complex128 *A, uni10_uint64 N);

    void setIdentity(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    // Lapack
    //
    void matrixSVD(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* U, uni10_complex128* S, uni10_complex128* vT);

    void matrixSDD(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* U, uni10_complex128* S, uni10_complex128* vT);

    void matrixQR(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* R);

    void matrixRQ(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* R, uni10_complex128* Q);

    void matrixQL(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* L);

    void matrixLQ(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* L, uni10_complex128* Q);

    void matrixQDR(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* D, uni10_complex128* R);

    void matrixLDQ(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* L, uni10_complex128* D, uni10_complex128* Q);

    void matrixQDRCPIVOT(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* Q, uni10_complex128* D, uni10_complex128* R);

    void matrixInv(uni10_complex128* A, uni10_int N);

    uni10_complex128 matrixDet(uni10_complex128* A, uni10_int N);

    //=================================================================================//

    void eigDecompose(uni10_complex128* Kij, uni10_int N, uni10_complex128* Eig, uni10_complex128 *EigVec);

    //
    // function overload for operators + - * += -= *=
    // 
    void matrix_diag_dense_add(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_sub(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_mul(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_add(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n);
                                                                                      
    void matrix_dense_diag_sub(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n);

    void matrix_diag_dense_mul(std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n);
   
    void matrix_dense_diag_add(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);
                                                                                                                                                    
    void matrix_dense_diag_sub(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);
                                                                                                                                                    
    void matrix_dense_diag_mul(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
