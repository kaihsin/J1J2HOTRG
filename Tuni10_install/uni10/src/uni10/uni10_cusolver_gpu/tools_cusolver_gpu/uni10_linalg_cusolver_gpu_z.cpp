#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_cusolver_gpu_z.h"

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd_gpu(std::complex<double> a, std::complex<double>* X, uni10_int incx, std::complex<double>* Y, uni10_int incy, uni10_uint64 N){   // Y = aX + Y
      
      uni10_error_msg(true, "%s", "Developing");

    }

    // Blas
    void vectorAdd_gpu(std::complex<double>* Y, std::complex<double>* X, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void vectorSub_gpu(std::complex<double>* Y, std::complex<double>* X, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void vectorMul_gpu(std::complex<double>* Y, std::complex<double>* X, uni10_uint64 N){ 

      uni10_error_msg(true, "%s", "Developing");

    }

    void vectorScal_gpu(std::complex<double> a, std::complex<double>* X, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void vectorExp_gpu(std::complex<double> a, std::complex<double>* X, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    std::complex<double> vectorSum_gpu(std::complex<double>* X, uni10_uint64 N, uni10_int inc){

      uni10_error_msg(true, "%s", "Developing");

    }

    double vectorNorm_gpu(std::complex<double>* X, uni10_uint64 N, uni10_int inc){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixDot_gpu(std::complex<double>* A, std::complex<double>* B, uni10_int M, uni10_int N, uni10_int K, std::complex<double>* C){

      uni10_error_msg(true, "%s", "Developing");

    }

    void diagRowMul_gpu(std::complex<double>* mat, std::complex<double>* diag, uni10_uint64 M, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void diagColMul_gpu(std::complex<double> *mat, std::complex<double>* diag, uni10_uint64 M, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void setTranspose_gpu(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N, std::complex<double>* AT){

      uni10_error_msg(true, "%s", "Developing");

    }

    void setTranspose_gpu(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void setDagger_gpu(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N, std::complex<double> *AT){

      uni10_error_msg(true, "%s", "Developing");

    }

    void setDagger_gpu(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void setConjugate_gpu(std::complex<double> *A, uni10_uint64 N, std::complex<double> *A_conj){

      uni10_error_msg(true, "%s", "Developing");

    }

    void setConjugate_gpu(std::complex<double> *A, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    //LAPACK
    //
    void matrixSVD_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, std::complex<double>* S_ori, std::complex<double>* vT){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixSDD_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, std::complex<double>* S_ori, std::complex<double>* vT){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixQR_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* R){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixRQ_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* R, std::complex<double>* Q){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixLQ_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* L, std::complex<double>* Q){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixQL_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* L){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixQDR_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* D, std::complex<double>* R){

      uni10_error_msg(true, "%s", "Developing");

    }


    void matrixLDQ_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* L, std::complex<double>* D, std::complex<double>* Q){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixQDRCPIVOT_gpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* D, std::complex<double>* R){

      uni10_error_msg(true, "%s", "Developing");

    }


    void matrixInv_gpu(std::complex<double>* A, uni10_int N){

      uni10_error_msg(true, "%s", "Developing");
    }


    std::complex<double> matrixDet_gpu(std::complex<double>* A, uni10_int N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void eigDecompose_gpu(std::complex<double>* Kij, uni10_int N, std::complex<double>* Eig, std::complex<double>* EigVec){
      uni10_uint64 memsize = N * N * sizeof(std::complex<double>);

      uni10_error_msg(true, "%s", "Developing");

    }

    void setIdentity_gpu(std::complex<double>* elem, uni10_uint64 M, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    //
    // function overload for operators + - * += -= *=
    // 
    void matrix_diag_dense_add_gpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_sub_gpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_mul_gpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_add_gpu(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_sub_gpu(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_mul_gpu(std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n){

      uni10_error_msg(true, "%s", "Developing");

    }
   
    void matrix_dense_diag_add_gpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_sub_gpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_mul_gpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_error_msg(true, "%s", "Developing");

    }

  } /* namespace uni10_linalg */

}
