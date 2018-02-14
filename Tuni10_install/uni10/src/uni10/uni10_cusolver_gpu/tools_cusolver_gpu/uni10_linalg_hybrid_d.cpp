#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_tools_cusolver_gpu.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_hybrid_d.h"

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(uni10_double64 a, uni10_double64* X, uni10_int incx, uni10_double64* Y, uni10_int incy, uni10_uint64 N){   // Y = aX + Y

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

    void vectorAdd(uni10_double64* Y, uni10_double64* X, uni10_uint64 N){   // Y = Y + X

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

    void vectorSub(uni10_double64* Y, uni10_double64* X, uni10_uint64 N){   // Y = Y + X

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

    void vectorMul(uni10_double64* Y, uni10_double64* X, uni10_uint64 N){ // Y = Y * X, element-wise multiplication;

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

    void vectorScal(double a, double* X, uni10_uint64 N){

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

    void vectorExp(double a, double* X, uni10_uint64 N){

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

    double vectorSum(double* X, uni10_uint64 N, uni10_int inc){

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

    double vectorNorm(double* X, uni10_uint64 N, uni10_int inc){

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

    void matrixDot(double* A, double* B, uni10_int M, uni10_int N, uni10_int K, double* C){

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

    void diagRowMul(double* mat, double* diag, uni10_uint64 M, uni10_uint64 N){

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

    void diagColMul(double *mat, double* diag, uni10_uint64 M, uni10_uint64 N){

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

    void setTranspose(double* A, uni10_uint64 M, uni10_uint64 N, double* AT){

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

    void setTranspose(double* A, uni10_uint64 M, uni10_uint64 N){

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

    void setDagger(double* A, uni10_uint64 M, uni10_uint64 N, double* AT){

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

    void setDagger(double* A, uni10_uint64 M, uni10_uint64 N){

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

    void matrixSVD(double* Mij_ori, uni10_int M, uni10_int N, double* U, double* S, double* vT){

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

    void matrixSDD(double* Mij_ori, uni10_int M, uni10_int N, double* U, double* S, double* vT){

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

    void matrixQR(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* R){

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

    void matrixRQ(double* Mij_ori, uni10_int M, uni10_int N, double* R, double* Q){

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

    void matrixLQ(double* Mij_ori, uni10_int M, uni10_int N, double* L, double* Q){

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

    void matrixQL(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* L){

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

    void matrixQDR(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* D, double* R){

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

    void matrixQDRCPIVOT(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* D, double* R){

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

    void matrixLDQ(double* Mij_ori, uni10_int M, uni10_int N, double* L, double* D, double* Q){

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

    void matrixInv(double* A, uni10_int N){

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

    double matrixDet(double* A, uni10_int N){

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

    void setIdentity(double* elem, uni10_uint64 M, uni10_uint64 N){

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

    void trimatrixEigH(double* D, double* E, uni10_int N, double* z, uni10_int LDZ){

      if(RUNTIMETYPE == only_cpu){

        trimatrixEigH_cpu(D, E, N, z, LDZ);

      }

      else if(RUNTIMETYPE == only_gpu){

        trimatrixEigH_gpu(D, E, N, z, LDZ);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }


    }

    void eigSyDecompose(double* Kij, uni10_int N, double* Eig, double* EigVec){

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

    //
    // function overload for operator + - * += -= *= 
    //
    void matrix_diag_dense_add(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b){

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

    void matrix_diag_dense_sub(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b){

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

    void matrix_diag_dense_mul(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b){

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

    void matrix_dense_diag_add(double* a, const double* D, uni10_uint64 m, uni10_uint64 n){

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

    void matrix_dense_diag_sub(double* a, const double* D, uni10_uint64 m, uni10_uint64 n){

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

    void matrix_diag_dense_mul(double* D, const double* a, uni10_uint64 m, uni10_uint64 n){

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

    void matrix_dense_diag_add(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b){

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

    void matrix_dense_diag_sub(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b){

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

    void matrix_dense_diag_mul(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b){

      if(RUNTIMETYPE == only_cpu){

        matrix_dense_diag_mul_cpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == only_gpu){

        matrix_dense_diag_mul_cpu(a, D, m, n, b);

      }

      else if(RUNTIMETYPE == hybrid){

        uni10_error_msg(true, "%s", "Developing");

      }

    }


  } /* namespace uni10_linalg */

} /* namespace uni10 */

