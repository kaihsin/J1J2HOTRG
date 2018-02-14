#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_cusolver_gpu_dz.h"

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd_gpu(double a, double* X, uni10_int incx, std::complex<double>* Y, uni10_int incy, uni10_uint64 N){   // Y = aX + Y

      uni10_error_msg(true, "%s", "Developing");

    }

    void vectorAdd_gpu(std::complex<double>* Y, double* X, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void vectorSub_gpu(std::complex<double>* Y, double* X, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void vectorMul_gpu(std::complex<double>* Y, double* X, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void vectorScal_gpu(double a, std::complex<double>* X, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void vectorExp_gpu(double a, std::complex<double>* X, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixDot_gpu(double* A, std::complex<double>* B, uni10_int M, uni10_int N, uni10_int K, std::complex<double>* C){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixDot_gpu(std::complex<double>* A, double* B, uni10_int M, uni10_int N, uni10_int K, std::complex<double>* C){

      uni10_error_msg(true, "%s", "Developing");

    }

    void diagRowMul_gpu(std::complex<double>* mat, double* diag, uni10_uint64 M, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void diagColMul_gpu(std::complex<double>* mat, double* diag, uni10_uint64 M, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void eigDecompose_gpu(double* Kij_ori, uni10_int N, std::complex<double>* Eig, std::complex<double>* EigVec){
      std::complex<double> *Kij = (std::complex<double>*) malloc(N * N * sizeof(std::complex<double>));

      uni10_error_msg(true, "%s", "Developing");

    }

    void eigSyDecompose_gpu(std::complex<double>* Kij, uni10_int N, double* Eig, std::complex<double>* EigVec){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixSVD_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, double *S, std::complex<double>* vT){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixSDD_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, double *S, std::complex<double>* vT){

      uni10_error_msg(true, "%s", "Developing");

    }

    //
    // function overload for operators + - * += -= *=
    // 
    // r operator(+,-,*,+=,-=,*=) z
    void matrix_diag_dense_add_gpu(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_sub_gpu(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_mul_gpu(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }
   
    void matrix_dense_diag_add_gpu(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_sub_gpu(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_mul_gpu(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    // z operator(+,-,*,+=,-=,*=) r

    void matrix_diag_dense_add_gpu(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_sub_gpu(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_mul_gpu(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_add_gpu(std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_sub_gpu(std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_mul_gpu(std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n){

      uni10_error_msg(true, "%s", "Developing");

    }
   
    void matrix_dense_diag_add_gpu(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_sub_gpu(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_mul_gpu(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }


  };/* namespace uni10_linalg */

};/* namespace uni10 */

