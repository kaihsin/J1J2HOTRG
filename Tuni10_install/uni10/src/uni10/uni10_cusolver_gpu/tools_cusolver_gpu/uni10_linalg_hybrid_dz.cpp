#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_tools_cusolver_gpu.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_hybrid_z.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_hybrid_dz.h"

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(double a, double* X, uni10_int incx, std::complex<double>* Y, uni10_int incy, uni10_uint64 N){   // Y = aX + Y

      if(RUNTIMETYPE == only_cpu){

        vectorAdd_cpu(a, X, incx, Y, incy, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        vectorAdd_gpu(a, X, incx, Y, incy, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void vectorAdd(std::complex<double>* Y, double* X, uni10_uint64 N){

      if(RUNTIMETYPE == only_cpu){

        vectorAdd_cpu(Y, X, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        vectorAdd_gpu(Y, X, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void vectorSub(std::complex<double>* Y, double* X, uni10_uint64 N){

      if(RUNTIMETYPE == only_cpu){

        vectorSub_cpu(Y, X, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        vectorSub_gpu(Y, X, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void vectorMul(std::complex<double>* Y, double* X, uni10_uint64 N){

      if(RUNTIMETYPE == only_cpu){

        vectorMul_cpu(Y, X, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        vectorMul_gpu(Y, X, N);


      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void vectorScal(double a, std::complex<double>* X, uni10_uint64 N){
     
      if(RUNTIMETYPE == only_cpu){

        vectorScal_cpu(a, X, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        vectorScal_gpu(a, X, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void vectorExp(double a, std::complex<double>* X, uni10_uint64 N){

      if(RUNTIMETYPE == only_cpu){

        vectorExp_cpu(a, X, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        vectorExp_gpu(a, X, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrixDot(double* A, std::complex<double>* B, uni10_int M, uni10_int N, uni10_int K, std::complex<double>* C){

      if(RUNTIMETYPE == only_cpu){

        matrixDot_cpu(A, B, M, N, K, C);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixDot_gpu(A, B, M, N, K, C);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrixDot(std::complex<double>* A, double* B, uni10_int M, uni10_int N, uni10_int K, std::complex<double>* C){

      if(RUNTIMETYPE == only_cpu){

        matrixDot_cpu(A, B, M, N, K, C);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixDot_gpu(A, B, M, N, K, C);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void diagRowMul(std::complex<double>* mat, double* diag, uni10_uint64 M, uni10_uint64 N){

      if(RUNTIMETYPE == only_cpu){

        diagRowMul_cpu(mat, diag, M, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        diagRowMul_gpu(mat, diag, M, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void diagColMul(std::complex<double>* mat, double* diag, uni10_uint64 M, uni10_uint64 N){

      if(RUNTIMETYPE == only_cpu){

        diagColMul_cpu(mat, diag, M, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        diagColMul_gpu(mat, diag, M, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void eigDecompose(double* Kij_ori, uni10_int N, std::complex<double>* Eig, std::complex<double>* EigVec){

      if(RUNTIMETYPE == only_cpu){

        eigDecompose_cpu(Kij_ori, N, Eig, EigVec);

      }

      else if(RUNTIMETYPE == only_gpu){

        eigDecompose_gpu(Kij_ori, N, Eig, EigVec);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void eigSyDecompose(std::complex<double>* Kij, uni10_int N, double* Eig, std::complex<double>* EigVec){

      if(RUNTIMETYPE == only_cpu){

        eigSyDecompose_cpu(Kij, N, Eig, EigVec);

      }

      else if(RUNTIMETYPE == only_gpu){

        eigSyDecompose_gpu(Kij, N, Eig, EigVec);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrixSVD(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, double *S, std::complex<double>* vT){

      if(RUNTIMETYPE == only_cpu){

        matrixSVD_cpu(Mij_ori, M, N, U, S, vT);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixSVD_gpu(Mij_ori, M, N, U, S, vT);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrixSDD(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, double *S, std::complex<double>* vT){

      if(RUNTIMETYPE == only_cpu){

        matrixSDD_cpu(Mij_ori, M, N, U, S, vT);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixSDD_gpu(Mij_ori, M, N, U, S, vT);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    //
    // function overload for operators + - * += -= *=
    // 
    // r operator(+,-,*,+=,-=,*=) z
    void matrix_diag_dense_add(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_diag_dense_add_cpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_diag_dense_add_gpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrix_diag_dense_sub(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_diag_dense_sub_cpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_diag_dense_sub_gpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrix_diag_dense_mul(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_diag_dense_mul_cpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_diag_dense_mul_gpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }
   
    void matrix_dense_diag_add(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_add_cpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_add_gpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }
                                                                                                                                      
    void matrix_dense_diag_sub(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_sub_cpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_sub_gpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }
                                                                                                                                      
    void matrix_dense_diag_mul(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_mul_cpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_mul_gpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    // z operator(+,-,*,+=,-=,*=) r

    void matrix_diag_dense_add(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_diag_dense_add_cpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_diag_dense_add_gpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrix_diag_dense_sub(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_sub_cpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_sub_gpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrix_diag_dense_mul(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_sub_cpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_sub_gpu(D, a, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrix_dense_diag_add(std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_add_cpu(a, D, m, n);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_add_gpu(a, D, m, n);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }
                                                                                                       
    void matrix_dense_diag_sub(std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_sub_cpu(a, D, m, n);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_sub_gpu(a, D, m, n);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }
                                                                                                       
    void matrix_diag_dense_mul(std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n){

      if(RUNTIMETYPE == only_cpu){

        matrix_diag_dense_mul_cpu(D, a, m, n);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_diag_dense_mul_gpu(D, a, m, n);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }
   
    void matrix_dense_diag_add(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_add_cpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_add_gpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrix_dense_diag_sub(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_sub_cpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_sub_gpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrix_dense_diag_mul(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_mul_cpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_mul_gpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }


  };/* namespace uni10_linalg */

};/* namespace uni10 */

