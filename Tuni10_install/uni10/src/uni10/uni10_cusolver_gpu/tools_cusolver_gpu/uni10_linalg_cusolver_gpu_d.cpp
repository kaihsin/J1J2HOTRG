#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_cusolver_gpu_d.h"

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd_gpu(uni10_double64 a, uni10_double64* X, uni10_int incx, uni10_double64* Y, uni10_int incy, uni10_uint64 N){   // Y = aX + Y

      int64_t left        = N;
      uni10_uint64 offset = 0;
      uni10_int chunk;
      while(left > 0){

        cublasHandle_t handle;
        checkCudaErrors(cublasCreate(&handle));

        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        checkCudaErrors(cublasDaxpy(handle, chunk, &a, X + offset, incx, Y + offset, incy));
        offset += chunk;
        left -= INT_MAX;

        checkCudaErrors(cublasDestroy(handle));

      }

    }

    void vectorAdd_gpu(uni10_double64* Y, uni10_double64* X, uni10_uint64 N){   // Y = Y + X

      uni10_double64 a    = 1.0;

      int64_t left        = N;
      uni10_uint64 offset = 0;
      uni10_int inc       = 1;
      uni10_int chunk;

      while(left > 0){

        cublasHandle_t handle;
        checkCudaErrors(cublasCreate(&handle));

        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        checkCudaErrors(cublasDaxpy(handle, chunk, &a, X + offset, inc, Y + offset, inc));
        offset += chunk;
        left -= INT_MAX;

        checkCudaErrors(cublasDestroy(handle));

      }


    }

    void vectorSub_gpu(uni10_double64* Y, uni10_double64* X, uni10_uint64 N){   // Y = Y + X


      uni10_double64 a    = -1.0;

      int64_t left        = N;
      uni10_uint64 offset = 0;
      uni10_int inc       = 1;
      uni10_int chunk;

      while(left > 0){

        cublasHandle_t handle;
        checkCudaErrors(cublasCreate(&handle));

        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        checkCudaErrors(cublasDaxpy(handle, chunk, &a, X + offset, inc, Y + offset, inc));
        offset += chunk;
        left -= INT_MAX;

        checkCudaErrors(cublasDestroy(handle));

      }

    }

    void vectorMul_gpu(double* Y, double* X, uni10_uint64 N){ // Y = Y * X, element-wise multiplication;

      vectorMul_kernel(Y, X, N);

    }

    void vectorScal_gpu(double a, double* X, uni10_uint64 N){

      int64_t left        = N;
      uni10_uint64 offset = 0;
      uni10_int inc       = 1;
      uni10_int chunk;

      while(left > 0){

        cublasHandle_t handle;
        checkCudaErrors(cublasCreate(&handle));

        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        checkCudaErrors(cublasDscal(handle, chunk, &a, X + offset, inc));
        offset += chunk;
        left -= INT_MAX;

        checkCudaErrors(cublasDestroy(handle));

      }

    }

    void vectorExp_gpu(double a, double* X, uni10_uint64 N){

      vectorExp_kernel(a, X, N);

    }

    double vectorSum_gpu(double* X, uni10_uint64 N, uni10_int inc){

      uni10_error_msg(true, "%s", "Developing");

    }

    double vectorNorm_gpu(double* X, uni10_uint64 N, uni10_int inc){

      double norm2 = 0;
      double tmp = 0;
      int64_t left = N;
      uni10_uint64 offset = 0;
      uni10_int chunk;

      while(left > 0){

        cublasHandle_t handle;
        checkCudaErrors(cublasCreate(&handle));

        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        checkCudaErrors(cublasDnrm2(handle, chunk, X + offset, inc, &tmp));
        norm2 += tmp * tmp;
        offset += chunk;
        left -= INT_MAX;

        checkCudaErrors(cublasDestroy(handle));

      }

      return sqrt(norm2);


    }

    void matrixDot_gpu(double* A, double* B, uni10_int M, uni10_int N, uni10_int K, double* C){

      double alpha = 1, beta = 0;

      cublasHandle_t handle;
      checkCudaErrors(cublasCreate(&handle));

      checkCudaErrors(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, N, M, K, &alpha, B, N, A, K, &beta, C, N));

      checkCudaErrors(cublasDestroy(handle));


    }

    void diagRowMul_gpu(double* mat, double* diag, uni10_uint64 M, uni10_uint64 N){

      double alpha = 1, beta = 0;
      uni10_int inc = 1;

      double* buffer;
      uni10_uint64 memsize = M * N * sizeof(double);

      checkCudaErrors(cudaMallocManaged(&buffer, memsize));

      cublasHandle_t handle;
      checkCudaErrors(cublasCreate(&handle));

      checkCudaErrors(cublasDdgmm(handle, CUBLAS_SIDE_RIGHT, N, M, mat, N, diag, inc, buffer, N));

      checkCudaErrors(cublasDestroy(handle));

      checkCudaErrors(cudaMemcpy(mat, buffer, memsize, cudaMemcpyDeviceToDevice));

      checkCudaErrors(cudaFree(buffer));

    }

    void diagColMul_gpu(double *mat, double* diag, uni10_uint64 M, uni10_uint64 N){

      double alpha = 1, beta = 0;
      uni10_int inc = 1;

      double* buffer;
      uni10_uint64 memsize = M * N * sizeof(double);

      checkCudaErrors(cudaMallocManaged(&buffer, memsize));

      cublasHandle_t handle;
      checkCudaErrors(cublasCreate(&handle));

      checkCudaErrors(cublasDdgmm(handle, CUBLAS_SIDE_LEFT, N, M, mat, N, diag, inc, buffer, N));

      checkCudaErrors(cublasDestroy(handle));

      checkCudaErrors(cudaMemcpy(mat, buffer, memsize, cudaMemcpyDeviceToDevice));

      checkCudaErrors(cudaFree(buffer));


    }

    void setTranspose_gpu(double* A, uni10_uint64 M, uni10_uint64 N, double* AT){

      double alpha = 1.0, beta = 0.0;
      double *dummyPtr = NULL;
      cublasHandle_t handle;
      checkCudaErrors(cublasCreate(&handle));

      checkCudaErrors(cublasDgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, M, N, &alpha, A, N, &beta, dummyPtr, M, AT, M));

      checkCudaErrors(cublasDestroy(handle));

    }

    void setTranspose_gpu(double* A, uni10_uint64 M, uni10_uint64 N){

      double* buffer;

      uni10_uint64 memsize = M * N * sizeof(double);

      checkCudaErrors(cudaMalloc(&buffer, memsize));
      checkCudaErrors(cudaMemcpy(buffer, A, memsize, cudaMemcpyDeviceToDevice));

      double alpha = 1.0, beta = 0.0;
      double *dummyPtr = NULL;

      cublasHandle_t handle;
      checkCudaErrors(cublasCreate(&handle));

      checkCudaErrors(cublasDgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, M, N, &alpha, buffer, N, &beta, dummyPtr, M, A, M));

      checkCudaErrors(cudaFree(buffer));
      checkCudaErrors(cublasDestroy(handle));

    }

    void setDagger_gpu(double* A, uni10_uint64 M, uni10_uint64 N, double* AT){

      setTranspose_gpu(A, M, N, AT);
  
    }

    void setDagger_gpu(double* A, uni10_uint64 M, uni10_uint64 N){

      setTranspose_gpu(A, M, N);

    }

    void matrixSVD_gpu(double* Mij_ori, uni10_int M, uni10_int N, double* U, double* S, double* vT){

      char jobu[1], jobv[1]; 
      jobu[0] = ( U  == NULL ) ? 'N' : 'S';
      jobv[0] = ( vT == NULL ) ? 'N' : 'S';

      cusolverDnHandle_t handle;
      checkCudaErrors(cusolverDnCreate(&handle));

      uni10_uint64 memsize = M * N * sizeof(double);

      double* Mij = NULL, *bufU = NULL, *bufvT = NULL;
      checkCudaErrors(cudaMalloc(&Mij, memsize));

      uni10_int op_m = N, op_n = M;
      uni10_int min = std::min(M, N);

      if(N < M){
        setTranspose_gpu(Mij_ori, M, N, Mij);
        checkCudaErrors(cudaMalloc(&bufU , M * min * sizeof(double)));
        checkCudaErrors(cudaMalloc(&bufvT, min * N * sizeof(double)));
        op_m = M;
        op_n = N;
      }
      else{
        checkCudaErrors(cudaMemcpy(Mij, Mij_ori, memsize, cudaMemcpyDeviceToDevice));
        bufU  = vT;
        bufvT = U;
      }

      uni10_int lwork = 0;
      double* rwork = NULL;
      double* work = NULL;
      uni10_int ldA = op_m, ldu = op_m, ldvT = min;
      checkCudaErrors(cusolverDnDgesvd_bufferSize(handle, op_m, op_n, &lwork));

      checkCudaErrors(cudaMalloc(&work,  sizeof(double)*lwork));
      checkCudaErrors(cudaMalloc(&rwork, sizeof(double)*(min-1)));

      uni10_int *devInfo;
      checkCudaErrors(cudaMalloc(&devInfo, sizeof(uni10_int)));
      checkCudaErrors(cudaMemset(devInfo, 0, sizeof(uni10_int)));

      checkCudaErrors(cusolverDnDgesvd(handle, jobv[0], jobu[0], op_m, op_n, Mij, ldA, S, bufU, ldu, bufvT, ldvT, work, lwork, rwork, devInfo));

      if(N < M){
        setTranspose_gpu(bufU,  min, M, U);
        setTranspose_gpu(bufvT, N, min, vT);
        checkCudaErrors(cudaFree(bufU));
        checkCudaErrors(cudaFree(bufvT));
      }

      uni10_int info;
      checkCudaErrors(cudaMemcpy(&info, devInfo, sizeof(uni10_int), cudaMemcpyDeviceToHost));
      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'cusolverDnDgesvd': cusolver DEVINFO = ", info);

      checkCudaErrors(cudaFree(Mij));
      checkCudaErrors(cudaFree(work));
      checkCudaErrors(cudaFree(rwork));
      checkCudaErrors(cudaFree(devInfo));
      checkCudaErrors(cusolverDnDestroy(handle));

    }

    void matrixSDD_gpu(double* Mij_ori, uni10_int M, uni10_int N, double* U, double* S, double* vT){

      uni10_error_msg(true, "%s", "Dose not support in cusolverDn-8.0. Try to use svd() instead.");

    }

    void matrixQR_gpu(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* R){

      double *Mij, *QT, *RT;
      checkCudaErrors(cudaMalloc(&Mij, M * N * sizeof(double)));
      checkCudaErrors(cudaMalloc(&QT , M * N * sizeof(double)));
      checkCudaErrors(cudaMalloc(&RT , N * N * sizeof(double)));

      setTranspose_gpu(Mij_ori, M, N, Mij);

      matrixLQ_gpu(Mij, N, M, RT, QT);

      setTranspose_gpu(QT, N, M, Q);
      setTranspose_gpu(RT, N, N, R);

      checkCudaErrors(cudaFree(Mij));
      checkCudaErrors(cudaFree(QT));
      checkCudaErrors(cudaFree(RT));

    }

    void matrixRQ_gpu(double* Mij_ori, uni10_int M, uni10_int N, double* R, double* Q){

      uni10_error_msg(true, "%s", "Dose not support in cusolverDn-8.0. Try to use qr() or lq() instead.");

    }

    void matrixLQ_gpu(double* Mij_ori, uni10_int M, uni10_int N, double* L, double* Q){

      uni10_uint64 memsize = M * N * sizeof(double);

      checkCudaErrors(cudaMemcpy(Q, Mij_ori, memsize, cudaMemcpyDeviceToDevice));
      
      double *tau = NULL, *work = NULL;
      uni10_int *devInfo = NULL;
      uni10_int lwork, lwork_geqrf, lwork_orgqr;

      uni10_int min = std::min(M, N);
      checkCudaErrors(cudaMalloc(&tau    , min * sizeof(double)));
      checkCudaErrors(cudaMalloc(&devInfo, sizeof(double)));

      cusolverDnHandle_t cusolverhandle;
      cusolverDnCreate(&cusolverhandle);

      cublasHandle_t cublashandle;
      cublasCreate(&cublashandle);

      checkCudaErrors(cusolverDnDgeqrf_bufferSize(cusolverhandle, N, M, Q, N, &lwork_geqrf));
      checkCudaErrors(cusolverDnDorgqr_bufferSize(cusolverhandle, N, M, M, Q, N, tau, &lwork_orgqr));

      lwork = (lwork_geqrf > lwork_orgqr) ? lwork_geqrf : lwork_orgqr;
      checkCudaErrors(cudaMalloc(&work   , lwork * sizeof(double)));

      checkCudaErrors(cusolverDnDgeqrf(cusolverhandle, N, M, Q, N, tau, work, lwork, devInfo));

      uni10_int info;
      checkCudaErrors(cudaMemcpy(&info, devInfo, sizeof(uni10_int), cudaMemcpyDeviceToHost));
      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'cusolverDnDgeqrf': cusolver DEVINFO = ", info);

      checkCudaErrors(cusolverDnDorgqr(cusolverhandle, N, M, M, Q, N, tau, work, lwork, devInfo));

      checkCudaErrors(cudaMemcpy(&info, devInfo, sizeof(uni10_int), cudaMemcpyDeviceToHost));
      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'cusolverDnDorgqr': cusolver DEVINFO = ", info);

      double alpha = 1, beta = 0;
      checkCudaErrors(cublasDgemm(cublashandle, CUBLAS_OP_T, CUBLAS_OP_N, M, M, N, &alpha, Q, N, Mij_ori, N, &beta, L, M));

      cublasDestroy(cublashandle);
      cusolverDnDestroy(cusolverhandle);

      checkCudaErrors(cudaFree(tau));
      checkCudaErrors(cudaFree(work));
      checkCudaErrors(cudaFree(devInfo));

    }

    void matrixQL_gpu(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* L){

      uni10_error_msg(true, "%s", "Dose not support in cusolverDn-8.0. Try to use qr() or lq() instead.");

    }

    void matrixQDR_gpu(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* D, double* R){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixQDRCPIVOT_gpu(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* D, double* R){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixLDQ_gpu(double* Mij_ori, uni10_int M, uni10_int N, double* L, double* D, double* Q){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrixInv_gpu(double* A, uni10_int N){

      uni10_error_msg(true, "%s", "Developing");

    }

    double matrixDet_gpu(double* A, uni10_int N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void setIdentity_gpu(double* elem, uni10_uint64 M, uni10_uint64 N){

      uni10_error_msg(true, "%s", "Developing");

    }

    void trimatrixEigH_gpu(double* D, double* E, uni10_int N, double* z, uni10_int LDZ){

      uni10_error_msg(true, "%s", "Developing");

    }

    void eigSyDecompose_gpu(double* Kij, uni10_int N, double* Eig, double* EigVec){

      checkCudaErrors(cudaMemcpy(EigVec, Kij, N * N * sizeof(double), cudaMemcpyDeviceToDevice));

      uni10_int lwork;
      cusolverDnHandle_t handle;
      
      checkCudaErrors(cusolverDnCreate(&handle));

      checkCudaErrors(cusolverDnDsyevd_bufferSize(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, N, EigVec, N, Eig, &lwork));

      double* work = NULL;
      uni10_int *devInfo = NULL;

      checkCudaErrors(cudaMalloc(&work   , lwork * sizeof(double)));
      checkCudaErrors(cudaMalloc(&devInfo, sizeof(uni10_int)));

      checkCudaErrors(cusolverDnDsyevd(handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, N, EigVec, N, Eig, work, lwork, devInfo));

      checkCudaErrors(cusolverDnDestroy(handle));

      checkCudaErrors(cudaFree(work));

    }

    //
    // function overload for operator + - * += -= *= 
    //
    void matrix_diag_dense_add_gpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_sub_gpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_mul_gpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_add_gpu(double* a, const double* D, uni10_uint64 m, uni10_uint64 n){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_sub_gpu(double* a, const double* D, uni10_uint64 m, uni10_uint64 n){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_diag_dense_mul_gpu(double* D, const double* a, uni10_uint64 m, uni10_uint64 n){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_add_gpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_sub_gpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b){

      uni10_error_msg(true, "%s", "Developing");

    }

    void matrix_dense_diag_mul_gpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b){

      uni10_error_msg(true, "%s", "Developing");

    }

  } /* namespace uni10_linalg */

} /* namespace uni10 */

