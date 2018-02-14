#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_lapack_cpu_z.h"


namespace uni10{

  namespace uni10_linalg{

    void vectorAdd_cpu(std::complex<double> a, std::complex<double>* X, uni10_int incx, std::complex<double>* Y, uni10_int incy, uni10_uint64 N){   // Y = aX + Y

      int64_t left      = N;
      uni10_uint64 offset = 0;
      uni10_int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        zaxpy(&chunk, &a, X + offset, &incx, Y + offset, &incy);
        offset += chunk;
        left -= INT_MAX;
      }
    }

    // Blas
    void vectorAdd_cpu(std::complex<double>* Y, std::complex<double>* X, uni10_uint64 N){

      std::complex<double> a = 1.0;

      int64_t left        = N;
      uni10_uint64 offset = 0;
      uni10_int inc       = 1;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        zaxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }

    }

    void vectorSub_cpu(std::complex<double>* Y, std::complex<double>* X, uni10_uint64 N){

      std::complex<double> a = -1.0;

      int64_t left      = N;
      uni10_uint64 offset = 0;
      uni10_int inc       = 1;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        zaxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }

    }

    void vectorMul_cpu(std::complex<double>* Y, std::complex<double>* X, uni10_uint64 N){ 
      for(uni10_uint64 i = 0; i < N; i++)
        Y[i] *= X[i];
    }

    void vectorScal_cpu(std::complex<double> a, std::complex<double>* X, uni10_uint64 N){

      int64_t left        = N;
      uni10_uint64 offset = 0;
      uni10_int inc       = 1;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        zscal(&chunk, &a, X + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }
    }

    void vectorExp_cpu(std::complex<double> a, std::complex<double>* X, uni10_uint64 N){
      for(uni10_uint64 i = 0; i < N; i++)
        X[i] = std::exp(a * X[i]);
    }

    std::complex<double> vectorSum_cpu(std::complex<double>* X, uni10_uint64 N, uni10_int inc){
      std::complex<double> sum = 0.0;
      uni10_uint64 idx = 0;
      for(uni10_uint64 i = 0; i < N; i++){
        sum += X[idx];
        idx += inc;
      }
      return sum;
    }

    double vectorNorm_cpu(std::complex<double>* X, uni10_uint64 N, uni10_int inc){

      double norm2 = 0;
      double tmp = 0;
      int64_t left = N;
      uni10_uint64 offset = 0;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        tmp = dznrm2(&chunk, X + offset, &inc);
        norm2 += tmp * tmp;
        offset += chunk;
        left -= INT_MAX;
      }
      return sqrt(norm2);

    }

    void matrixDot_cpu(std::complex<double>* A, std::complex<double>* B, uni10_int M, uni10_int N, uni10_int K, std::complex<double>* C){
      std::complex<double> alpha = 1.0, beta = 0.0;
      zgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);
    }

    void diagRowMul_cpu(std::complex<double>* mat, std::complex<double>* diag, uni10_uint64 M, uni10_uint64 N){
      for(uni10_uint64 i = 0; i < M; i++)
        vectorScal_cpu(diag[i], &(mat[i * N]), N);
    }

    void diagColMul_cpu(std::complex<double> *mat, std::complex<double>* diag, uni10_uint64 M, uni10_uint64 N){
      for(uni10_uint64 i = 0; i < M; i++){
        uni10_uint64 ridx = i * N;
        for(uni10_uint64 j = 0; j < N; j++)
          mat[ridx + j] *= diag[j];
      }
    }

    void setTranspose_cpu(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N, std::complex<double>* AT){
      for(uni10_uint64 i = 0; i < M; i++)
        for(uni10_uint64 j = 0; j < N; j++)
          AT[j * M + i] = A[i * N + j];
    }

    void setTranspose_cpu(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N){
      uni10_uint64 memsize = M * N * sizeof(std::complex<double>);
      std::complex<double> *AT = (std::complex<double>*)malloc(memsize);
      setTranspose_cpu(A, M, N, AT);
      memcpy(A, AT, memsize);
      free(AT);
    }

    void setDagger_cpu(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N, std::complex<double> *AT){
      for(uni10_uint64 i = 0; i < M; i++)
        for(uni10_uint64 j = 0; j < N; j++)
          AT[j * M + i] = std::conj(A[i * N + j]);
    }

    void setDagger_cpu(std::complex<double>* A, uni10_uint64 M, uni10_uint64 N){
      uni10_uint64 memsize = M * N * sizeof(std::complex<double>);
      std::complex<double> *AT = (std::complex<double>*)malloc(memsize);
      setDagger_cpu(A, M, N, AT);
      memcpy(A, AT, memsize);
      free(AT);
    }

    void setConjugate_cpu(std::complex<double> *A, uni10_uint64 N, std::complex<double> *A_conj){
      for(uni10_uint64 i = 0; i < N; i++)
        A_conj[i] = std::conj(A[i]);
    }

    void setConjugate_cpu(std::complex<double> *A, uni10_uint64 N){
      for(uni10_uint64 i = 0; i < N; i++)
        A[i] = std::conj(A[i]);
    }

    //LAPACK
    //
    void matrixSVD_cpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, std::complex<double>* S_ori, std::complex<double>* vT){

      uni10_int min = std::min(M, N);
      double* S = (double*)malloc(min * sizeof(double));
      matrixSVD_cpu(Mij_ori, M, N, U, S, vT);
      uni10_elem_cast(S_ori, S, min);
      free(S);

    }

    void matrixSDD_cpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, std::complex<double>* S_ori, std::complex<double>* vT){

      uni10_int min = std::min(M, N);
      double* S = (double*)malloc(min * sizeof(double));
      matrixSDD_cpu(Mij_ori, M, N, U, S, vT);
      uni10_elem_cast(S_ori, S, min);
      free(S);

    }

    void matrixQR_cpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* R){

      memcpy(Q, Mij_ori, N*M*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      uni10_int lda = N;
      uni10_int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      uni10_int info;
      uni10_int K = N;
      zgelqf(&N, &M, Q, &lda, tau, &worktestzge, &lwork, &info);
      zunglq(&N, &M, &K, Q, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (uni10_int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgelqf(&N, &M, Q, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getUpTri(Q, R, M, N);
      //getQ
      lwork = (uni10_int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zunglq(&N, &M, &K, Q, &lda, tau, workzun, &lwork, &info);
      
      //std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
      //zgemm((char*)"N", (char*)"C", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

      free(tau);
      free(workzge);
      free(workzun);
    }

    void matrixRQ_cpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* R, std::complex<double>* Q){

      memcpy(Q, Mij_ori, M*N*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      uni10_int lda = N;
      uni10_int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      uni10_int info;
      uni10_int K = M;
      zgeqlf(&N, &M, Q, &lda, tau, &worktestzge, &lwork, &info);
      zungql(&N, &M, &K, Q, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (uni10_int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgeqlf(&N, &M, Q, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getUpTri(Q, R, M, N);
      //getQ
      lwork = (uni10_int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zungql(&N, &M, &K, Q, &lda, tau, workzun, &lwork, &info);

      //std::complex<double> alpha (1.0, 0.0), beta (0.0, 0.0);
      //zgemm((char*)"C", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, R, &M);

      free(tau);
      free(workzge);
      free(workzun);

    }

    void matrixLQ_cpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* L, std::complex<double>* Q){

      memcpy(Q, Mij_ori, M*N*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      uni10_int lda = N;
      uni10_int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      uni10_int info;
      uni10_int K = M;
      zgeqrf(&N, &M, Q, &lda, tau, &worktestzge, &lwork, &info);
      zungqr(&N, &M, &K, Q, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (uni10_int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgeqrf(&N, &M, Q, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getDnTri(Q, L, M, N);
      //getQ
      lwork = (uni10_int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zungqr(&N, &M, &K, Q, &lda, tau, workzun, &lwork, &info);

      //std::complex<double> alpha (1.0, 0.0), beta (0.0, 0.0);
      //zgemm((char*)"C", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, L, &M);

      free(tau);
      free(workzge);
      free(workzun);

    }

    void matrixQL_cpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* L){

      memcpy(Q, Mij_ori, N*M*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      uni10_int lda = N;
      uni10_int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      uni10_int info;
      uni10_int K = N;
      zgerqf(&N, &M, Q, &lda, tau, &worktestzge, &lwork, &info);
      zungrq(&N, &M, &K, Q, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (uni10_int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgerqf(&N, &M, Q, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getDnTri(Q, L, M, N);
      //getQ
      lwork = (uni10_int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zungrq(&N, &M, &K, Q, &lda, tau, workzun, &lwork, &info);

      //std::complex<double> alpha (1.0, 0.0), beta (1.0, 1.0);
      //zgemm((char*)"N", (char*)"C", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, L, &N);

      free(tau);
      free(workzge);
      free(workzun);
    }

    void matrixQDR_cpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* D, std::complex<double>* R){

      memcpy(Q, Mij_ori, N*M*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      uni10_int lda = N;
      uni10_int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      uni10_int info;
      uni10_int K = N;
      zgelqf(&N, &M, Q, &lda, tau, &worktestzge, &lwork, &info);
      zunglq(&N, &M, &K, Q, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (uni10_int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgelqf(&N, &M, Q, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getUpTri(Q, R, M, N);
      uni10_getDiag(R, D, N, N, N);
      for(uni10_int i = 0; i < N; i++)
        for(uni10_int j = 0; j < N-i; j++)
          R[i*N+i+j] /= D[i];
      //getQ
      lwork = (uni10_int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zunglq(&N, &M, &K, Q, &lda, tau, workzun, &lwork, &info);

      //std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);
      //zgemm((char*)"N", (char*)"C", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

      free(tau);
      free(workzge);
      free(workzun);
    }


    void matrixLDQ_cpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* L, std::complex<double>* D, std::complex<double>* Q){

      memcpy(Q, Mij_ori, M*N*sizeof(std::complex<double>));
      std::complex<double>* tau = (std::complex<double>*)malloc(M*sizeof(std::complex<double>));
      uni10_int lda = N;
      uni10_int lwork = -1;
      std::complex<double> worktestzge;
      std::complex<double> worktestzun;
      uni10_int info;
      uni10_int K = M;
      zgeqrf(&N, &M, Q, &lda, tau, &worktestzge, &lwork, &info);
      zungqr(&N, &M, &K, Q, &lda, tau, &worktestzun, &lwork, &info);
      lwork = (uni10_int)worktestzge.real();
      std::complex<double>* workzge = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgeqrf(&N, &M, Q, &lda, tau, workzge, &lwork, &info);
      //getR
      uni10_getDnTri(Q, L, M, N);
      uni10_getDiag(L, D, M, M, M);
      for(uni10_int i = 0; i < M; i++)
        for(uni10_int j = 0; j < M-i; j++)
          L[(i+j)*M+i] /= D[i];
      //getQ
      lwork = (uni10_int)worktestzun.real();
      std::complex<double>* workzun = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zungqr(&N, &M, &K, Q, &lda, tau, workzun, &lwork, &info);

      //std::complex<double> alpha (1.0, 0.0), beta (0.0, 0.0);
      //zgemm((char*)"C", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, L, &M);

      free(tau);
      free(workzge);
      free(workzun);

    }

    void matrixQDRCPIVOT_cpu(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* Q, std::complex<double>* D, std::complex<double>* R){

      uni10_error_msg(M != N, "%s", "M must be equalt to N");

      std::complex<double>* Mij = (std::complex<double>*)malloc(M * N * sizeof(std::complex<double>)); // Q(M x M): orthogonal basis
      setTranspose_cpu(Mij_ori, M, N, Mij);                      //column major + square matrix  // TP(M x M)
      uni10_int max = M > N ? M : N;
      uni10_int min = M < N ? M : N;
      uni10_int lda = max;
      uni10_int lrwork = 2* N;
      double* rwork = (double*) malloc(lrwork * sizeof(double));
      uni10_int* jpvt = (uni10_int*)malloc(N * sizeof(uni10_int));  //column vectors
      memset(jpvt, 0, N * sizeof(uni10_int));
      std::complex<double>* tau = (std::complex<double>*)malloc(min * sizeof(std::complex<double>));
      uni10_int info;
      uni10_int lwork = -1;
      std::complex<double> worktest;
      zgeqp3(&M, &N, Mij, &lda, jpvt, tau, &worktest, &lwork, rwork, &info);
      lwork = (uni10_int)worktest.real();
      std::complex<double>* work = (std::complex<double>*) malloc(lwork * sizeof(std::complex<double>));
      zgeqp3(&M, &N, Mij, &lda, jpvt, tau, work, &lwork, rwork, &info);
      uni10_error_msg(info != 0, "Lapack Info = %d", info);
      for(uni10_int i = 0; i < M; i++)
        D[i] = Mij[i * N + i];                               // D
      std::complex<double>* T = (std::complex<double>*)calloc(M * N, sizeof(std::complex<double>));
      for(uni10_int i = 0; i < M; i++)
        for(uni10_int j = i; j < N; j++)
          if(i == j)
            T[i * N + j] = 1;
          else
            T[i * N + j] = Mij[j * N + i] / D[i];
      for(uni10_int i = 0; i < M; i++)
        for(uni10_int j = 0; j < N; j++)
          R[i * N + (jpvt[j]-1)] = T[i * N + j];              // R 
      zungqr(&M, &N, &N, Mij, &lda, tau, work, &lwork, &info);
      uni10_error_msg(info != 0, "Lapack Info = %d", info);
      setTranspose_cpu(Mij, M, N, Q);                             // Q
      free(Mij);
      free(rwork);
      free(work);
      free(jpvt);
      free(tau);
      free(T);

    }


    void matrixInv_cpu(std::complex<double>* A, uni10_int N){
      //if(diag){
      //  for(uni10_int i = 0; i < N; i++)
      //    A[i] = std::abs(A[i]) == 0 ? 0.0 : 1.0/A[i];
      //  return;
      //}
      uni10_int *ipiv = (uni10_int*)malloc((N+1) * sizeof(uni10_int));
      uni10_int info;
      zgetrf(&N, &N, A, &N, ipiv, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zgetrf': Lapack INFO = ", info);

      uni10_int lwork = -1;
      std::complex<double> worktest;
      zgetri(&N, A, &N, ipiv, &worktest, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zgetri': Lapack INFO = ", info);

      lwork = (uni10_int)(worktest.real());
      std::complex<double> *work = (std::complex<double>*)malloc(lwork * sizeof(std::complex<double>));
      zgetri(&N, A, &N, ipiv, work, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zgetri': Lapack INFO = ", info);

      free(ipiv);
      free(work);
    }

    std::complex<double> matrixDet_cpu(std::complex<double>* A, uni10_int N){

      uni10_int *ipiv = (uni10_int*)malloc((N+1)*sizeof(uni10_int));
      uni10_int lwork = 64 * N;
      std::complex<double> *work = (std::complex<double>*)malloc(lwork * sizeof(std::complex<double>));
      uni10_int info;
      zgetrf(&N,&N,A,&N,ipiv,&info);
      uni10_error_msg( info != 0, "%s %d", "Error in Lapack function 'zgetrf': Lapack INFO = ", info );
      std::complex<double> det = 1;
      uni10_int neg = 0;
      for (uni10_int i = 0; i < N; i++) {
        det *= A[i * N + i];
        if (ipiv[i] != (i+1)) neg = !neg;
      }
      free(ipiv);
      free(work);
      return neg?-det:det;

    }

    void eigDecompose_cpu(std::complex<double>* Kij, uni10_int N, std::complex<double>* Eig, std::complex<double>* EigVec){
      uni10_uint64 memsize = N * N * sizeof(std::complex<double>);
      std::complex<double> *A = (std::complex<double>*) malloc(memsize);
      memcpy(A, Kij, memsize);
      uni10_int ldA = N;
      uni10_int ldvl = 1;
      uni10_int ldvr = N;
      uni10_int lwork = -1;
      double *rwork = (double*) malloc(2 * N * sizeof(double));
      std::complex<double> worktest;
      uni10_int info;
      zgeev((char*)"N", (char*)"V", &N, A, &ldA, Eig, NULL, &ldvl, EigVec, &ldvr, &worktest, &lwork, rwork, &info);

      uni10_error_msg(info != 0, "%s, %d", "Error in Lapack function 'zgeev': Lapack INFO = ", info);

      lwork = (uni10_int)worktest.real();
      std::complex<double>* work = (std::complex<double>*)malloc(sizeof(std::complex<double>)*lwork);
      zgeev((char*)"N", (char*)"V", &N, A, &ldA, Eig, NULL, &ldvl, EigVec, &ldvr, work, &lwork, rwork, &info);

      uni10_error_msg(info != 0, "%s, %d", "Error in Lapack function 'zgeev': Lapack INFO = ", info);

      free(work);
      free(rwork);
      free(A);
    }

    void setIdentity_cpu(std::complex<double>* elem, uni10_uint64 M, uni10_uint64 N){
      uni10_uint64 min;
      if(M < N)
        min = M;
      else
        min = N;
      memset(elem, 0, M * N * sizeof(std::complex<double>));
      for(uni10_uint64 i = 0; i < min; i++)
        elem[i * N + i] = 1.0;
    }

    void matrix_diag_dense_add_cpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      memcpy(b, a, m*n*sizeof(uni10_complex128));
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }

    void matrix_diag_dense_sub_cpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_uint64 elemNum = m*n;
      memcpy(b, a, elemNum*sizeof(uni10_complex128));
      vectorScal_cpu(-1., b,elemNum);
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }

    void matrix_diag_dense_mul_cpu(const std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* v){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        v[i] = D[i] * a[i*n+i];

    }

    void matrix_dense_diag_add_cpu(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        a[i*n+i] += D[i];

    }
                                                                                      
    void matrix_dense_diag_sub_cpu(std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        a[i*n+i] -= D[i];

    }
                                                                                      
    void matrix_diag_dense_mul_cpu(std::complex<double>* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        D[i] *= a[i*n+i];

    }
   
    void matrix_dense_diag_add_cpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      memcpy(b, a, m*n*sizeof(uni10_complex128));
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }
                                                                                                                                                    
    void matrix_dense_diag_sub_cpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      memcpy(b, a, m*n*sizeof(uni10_complex128));
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] -= D[i];

    }
                                                                                                                                                    
    void matrix_dense_diag_mul_cpu(const std::complex<double>* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* v){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        v[i] = D[i] * a[i*n+i];

    }


  } /* namespace uni10_linalg */

}
