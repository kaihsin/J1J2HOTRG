#include <limits.h>

#include "uni10/uni10_error.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_tools_cusolver_gpu.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_hybrid_z.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_hybrid_dz.h"

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(std::complex<double> a, std::complex<double>* X, uni10_int incx, std::complex<double>* Y, uni10_int incy, uni10_uint64 N){   // Y = aX + Y

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

    // Blas
    void vectorAdd(std::complex<double>* Y, std::complex<double>* X, uni10_uint64 N){

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

    void vectorSub(std::complex<double>* Y, std::complex<double>* X, uni10_uint64 N){

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

    void vectorMul(std::complex<double>* Y, std::complex<double>* X, uni10_uint64 N){ 

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

    void vectorScal(std::complex<double> a, std::complex<double>* X, uni10_uint64 N){

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

    void vectorExp(std::complex<double> a, std::complex<double>* X, uni10_uint64 N){

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

    std::complex<double> vectorSum(std::complex<double>* X, uni10_uint64 N, uni10_int inc){

      if(RUNTIMETYPE == only_cpu){

        vectorSum_cpu(X, N, inc);

      }

      else if(RUNTIMETYPE == only_gpu){

        vectorSum_gpu(X, N, inc);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    double vectorNorm(std::complex<double>* X, uni10_uint64 N, uni10_int inc){

      if(RUNTIMETYPE == only_cpu){

        vectorNorm_cpu(X, N, inc);

      }

      else if(RUNTIMETYPE == only_gpu){

        vectorNorm_gpu(X, N, inc);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrixDot(std::complex<double>* A, std::complex<double>* B, uni10_int M, uni10_int N, uni10_int K, std::complex<double>* C){

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

    void diagRowMul(std::complex<double>* mat, std::complex<double>* diag, uni10_uint64 M, uni10_uint64 N){

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

    void diagColMul(std::complex<double> *mat, std::complex<double>* diag, uni10_uint64 M, uni10_uint64 N){

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

    void setTranspose(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N, std::complex<double>* AT){

      if(RUNTIMETYPE == only_cpu){

        setTranspose_cpu(A, M, N, AT);

      }

      else if(RUNTIMETYPE == only_gpu){

        setTranspose_gpu(A, M, N, AT);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void setTranspose(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N){

      if(RUNTIMETYPE == only_cpu){

        setTranspose_cpu(A, M, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        setTranspose_gpu(A, M, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void setDagger(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N, std::complex<double> *AT){

      if(RUNTIMETYPE == only_cpu){

        setDagger_cpu(A, M, N, AT);

      }

      else if(RUNTIMETYPE == only_gpu){

        setDagger_gpu(A, M, N, AT);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void setDagger(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N){

      if(RUNTIMETYPE == only_cpu){

        setDagger_cpu(A, M, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        setDagger_gpu(A, M, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void setConjugate(std::complex<double> *A, uni10_uint64 N, std::complex<double> *A_conj){

      if(RUNTIMETYPE == only_cpu){

        setConjugate_cpu(A, N, A_conj);

      }

      else if(RUNTIMETYPE == only_gpu){

        setConjugate_gpu(A, N, A_conj);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void setConjugate(std::complex<double> *A, uni10_uint64 N){

      if(RUNTIMETYPE == only_cpu){

        setConjugate_cpu(A, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        setConjugate_gpu(A, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    //LAPACK
    //
    void matrixSVD(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, std::complex<double>* S, std::complex<double>* vT){

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

    void matrixSDD(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, std::complex<double>* S, std::complex<double>* vT){

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

    void matrixQR(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* R){

      if(RUNTIMETYPE == only_cpu){

        matrixQR_cpu(Mij_ori, M, N, Q, R);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixQR_gpu(Mij_ori, M, N, Q, R);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrixRQ(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* R, std::complex<double>* Q){

      if(RUNTIMETYPE == only_cpu){

        matrixRQ_cpu(Mij_ori, M, N, R, Q);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixRQ_gpu(Mij_ori, M, N, R, Q);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrixLQ(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* L, std::complex<double>* Q){

      if(RUNTIMETYPE == only_cpu){

        matrixLQ_cpu(Mij_ori, M, N, L, Q);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixLQ_gpu(Mij_ori, M, N, L, Q);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrixQL(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* L){

      if(RUNTIMETYPE == only_cpu){

        matrixQL_cpu(Mij_ori, M, N, Q, L);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixQL_gpu(Mij_ori, M, N, Q, L);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrixQDR(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* D, std::complex<double>* R){

      if(RUNTIMETYPE == only_cpu){

        matrixQDR_cpu(Mij_ori, M, N, Q, D, R);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixQDR_gpu(Mij_ori, M, N, Q, D, R);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }


    void matrixLDQ(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* L, std::complex<double>* D, std::complex<double>* Q){

      if(RUNTIMETYPE == only_cpu){

        matrixLDQ_cpu(Mij_ori, M, N, L, D, Q);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixLDQ_gpu(Mij_ori, M, N, L, D, Q);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void matrixQDRCPIVOT(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* D, std::complex<double>* R){

      if(RUNTIMETYPE == only_cpu){

        matrixQDRCPIVOT_cpu(Mij_ori, M, N, Q, D, R);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixQDRCPIVOT_gpu(Mij_ori, M, N, Q, D, R);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }


    void matrixInv(std::complex<double>* A, uni10_int N){

      if(RUNTIMETYPE == only_cpu){

        matrixInv_cpu(A, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixInv_gpu(A, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    std::complex<double> matrixDet(std::complex<double>* A, uni10_int N){

      if(RUNTIMETYPE == only_cpu){

        matrixDet_cpu(A, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrixDet_gpu(A, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void eigDecompose(std::complex<double>* Kij, uni10_int N, std::complex<double>* Eig, std::complex<double>* EigVec){

      if(RUNTIMETYPE == only_cpu){

        eigDecompose_cpu(Kij, N, Eig, EigVec);

      }

      else if(RUNTIMETYPE == only_gpu){

        eigDecompose_gpu(Kij, N, Eig, EigVec);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    void setIdentity(std::complex<double>* elem, uni10_uint64 M, uni10_uint64 N){

      if(RUNTIMETYPE == only_cpu){

        setIdentity_cpu(elem, M, N);

      }

      else if(RUNTIMETYPE == only_gpu){

        setIdentity_gpu(elem, M, N);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }

    //
    // function overload for operators + - * += -= *=
    // 
    void matrix_diag_dense_add(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

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

    void matrix_diag_dense_sub(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

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

    void matrix_diag_dense_mul(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

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

    void matrix_dense_diag_add(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n){

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
                                                                                      
    void matrix_dense_diag_sub(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n){

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

    void matrix_diag_dense_mul(std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n){

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
   
    void matrix_dense_diag_add(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

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
                                                                                                                                                    
    void matrix_dense_diag_sub(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

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
                                                                                                                                                    
    void matrix_dense_diag_mul(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

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

  } /* namespace uni10_linalg */

}
