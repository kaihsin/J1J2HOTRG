#ifndef __UNI10_LINALG_LAPACK_CPU_DZ_H__
#define __UNI10_LINALG_LAPACK_CPU_DZ_H__

#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_tools_lapack_cpu.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_linalg_lapack_cpu_z.h"

#ifdef UNI_MKL
#include "mkl.h"
#else
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_lapack_wrapper_cpu.h"
#endif

namespace uni10{

  namespace uni10_linalg{

    //Blas
    //
    void vectorAdd(uni10_double64 a, uni10_double64* X, uni10_int incx, uni10_complex128* Y, uni10_int incy, uni10_uint64 N);   // Y = a*X + Y

    void vectorAdd(uni10_complex128* Y, uni10_double64* X, uni10_uint64 N);

    void vectorSub(uni10_complex128* Y, uni10_double64* X, uni10_uint64 N);

    void vectorMul(uni10_complex128* Y, uni10_double64* X, uni10_uint64 N);
    
    void vectorScal(uni10_double64 a, uni10_complex128* X, uni10_uint64 N);

    void vectorExp(uni10_double64 a, uni10_complex128* X, uni10_uint64 N);
 
    void matrixDot(uni10_double64* A, uni10_complex128* B, uni10_int M, uni10_int N, uni10_int K, uni10_complex128* C);

    void matrixDot(uni10_complex128* A, uni10_double64* B, uni10_int M, uni10_int N, uni10_int K, uni10_complex128* C);

    void diagRowMul(uni10_complex128* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul(uni10_complex128* mat, uni10_double64* diag, uni10_uint64 M, uni10_uint64 N);

    //Lapack
    //

    //====================================================================//

    void eigDecompose(uni10_double64* Kij, uni10_int N, uni10_complex128* Eig, uni10_complex128 *EigVec);

    void eigSyDecompose(uni10_complex128* Kij, uni10_int N, uni10_double64* Eig, uni10_complex128* EigVec);

    void matrixSVD(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* U, uni10_double64 *S, uni10_complex128* vT);

    void matrixSDD(uni10_complex128* Mij_ori, uni10_int M, uni10_int N, uni10_complex128* U, uni10_double64 *S, uni10_complex128* vT);
   
    //
    // function overload for operators + - * += -= *=
    // 
    // r operator(+,-,*,+=,-=,*=) z
    void matrix_diag_dense_add(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_sub(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_mul(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);
   
    void matrix_dense_diag_add(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);
                                                                                                                                      
    void matrix_dense_diag_sub(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);
                                                                                                                                      
    void matrix_dense_diag_mul(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    // z operator(+,-,*,+=,-=,*=) r

    void matrix_diag_dense_add(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_sub(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_diag_dense_mul(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

    void matrix_dense_diag_add(std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n);
                                                                                                       
    void matrix_dense_diag_sub(std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n);
                                                                                                       
    void matrix_diag_dense_mul(std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n);
   
    void matrix_dense_diag_add(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);
                                                                                                                                       
    void matrix_dense_diag_sub(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);
                                                                                                                                       
    void matrix_dense_diag_mul(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
