#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_linalg_lapack_cpu_d.h"


namespace uni10{

  namespace uni10_linalg{

    void vectorAdd_cpu(uni10_double64 a, uni10_double64* X, uni10_int incx, uni10_double64* Y, uni10_int incy, uni10_uint64 N){   // Y = aX + Y

      int64_t left        = N;
      uni10_uint64 offset = 0;
      uni10_int chunk;
      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        daxpy(&chunk, &a, X + offset, &incx, Y + offset, &incy);
        offset += chunk;
        left -= INT_MAX;
      }

    }

    void vectorAdd_cpu(uni10_double64* Y, uni10_double64* X, uni10_uint64 N){   // Y = Y + X

      uni10_double64 a    = 1.0;

      int64_t left        = N;
      uni10_uint64 offset = 0;
      uni10_int inc       = 1;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        daxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }

    }

    void vectorSub_cpu(uni10_double64* Y, uni10_double64* X, uni10_uint64 N){   // Y = Y + X

      uni10_double64 a    = -1.0;

      int64_t left        = N;
      uni10_uint64 offset = 0;
      uni10_int inc       = 1;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        daxpy(&chunk, &a, X + offset, &inc, Y + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }

    }

    void vectorMul_cpu(double* Y, double* X, uni10_uint64 N){ // Y = Y * X, element-wise multiplication;
      for(uni10_uint64 i = 0; i < N; i++)
        Y[i] *= X[i];
    }

    void vectorScal_cpu(double a, double* X, uni10_uint64 N){

      int64_t left        = N;
      uni10_uint64 offset = 0;
      uni10_int inc       = 1;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        dscal(&chunk, &a, X + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }
    }

    void vectorExp_cpu(double a, double* X, uni10_uint64 N){

      for(uni10_uint64 i = 0; i < N; i++)
        X[i] = std::exp(a * X[i]);

    }

    double vectorSum_cpu(double* X, uni10_uint64 N, uni10_int inc){

      double sum = 0;
      uni10_uint64 idx = 0;
      for(uni10_uint64 i = 0; i < N; i++){
        sum += X[idx];
        idx += inc;
      }
      return sum;

    }

    double vectorNorm_cpu(double* X, uni10_uint64 N, uni10_int inc){

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
        tmp = dnrm2(&chunk, X + offset, &inc);
        norm2 += tmp * tmp;
        offset += chunk;
        left -= INT_MAX;
      }

      return sqrt(norm2);

    }

    void matrixDot_cpu(double* A, double* B, uni10_int M, uni10_int N, uni10_int K, double* C){
      double alpha = 1, beta = 0;
      dgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);
    }

    void diagRowMul_cpu(double* mat, double* diag, uni10_uint64 M, uni10_uint64 N){

      for(uni10_uint64 i = 0; i < M; i++)
        vectorScal_cpu(diag[i], &(mat[i * N]), N);

    }

    void diagColMul_cpu(double *mat, double* diag, uni10_uint64 M, uni10_uint64 N){

      for(uni10_uint64 i = 0; i < M; i++){
        uni10_uint64 ridx = i * N;
        for(uni10_uint64 j = 0; j < N; j++)
          mat[ridx + j] *= diag[j];
      }

    }

    void setTranspose_cpu(double* A, uni10_uint64 M, uni10_uint64 N, double* AT){
      for(uni10_uint64 i = 0; i < M; i++)
        for(uni10_uint64 j = 0; j < N; j++)
          AT[j * M + i] = A[i * N + j];
    }

    void setTranspose_cpu(double* A, uni10_uint64 M, uni10_uint64 N){
      uni10_uint64 memsize = M * N * sizeof(double);
      double *AT = (double*)malloc(memsize);
      setTranspose_cpu(A, M, N, AT);
      memcpy(A, AT, memsize);
      free(AT);
    }

    void setDagger_cpu(double* A, uni10_uint64 M, uni10_uint64 N, double* AT){
      setTranspose_cpu(A, M, N, AT);
    }

    void setDagger_cpu(double* A, uni10_uint64 M, uni10_uint64 N){
      setTranspose_cpu(A, M, N);
    }

    void matrixSVD_cpu(double* Mij_ori, uni10_int M, uni10_int N, double* U, double* S, double* vT){

      char jobu[1], jobv[1]; 
      jobu[0] = ( U  == NULL ) ? 'N' : 'S';
      jobv[0] = ( vT == NULL ) ? 'N' : 'S';

      double* Mij = (double*)malloc(M * N * sizeof(double));
      memcpy(Mij, Mij_ori, M * N * sizeof(double));
      uni10_int min = std::min(M, N);
      uni10_int ldA = N, ldu = N, ldvT = min;
      uni10_int lwork = -1;
      double worktest;
      uni10_int info;

      dgesvd(jobv, jobu, &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, &worktest, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgesvd': Lapack INFO = ", info);

      lwork = (uni10_int)worktest;
      double *work = (double*)malloc(lwork*sizeof(double));
      dgesvd(jobv, jobu, &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, work, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgesvd': Lapack INFO = ", info);

      free(work);
      free(Mij);
    }

    void matrixSDD_cpu(double* Mij_ori, uni10_int M, uni10_int N, double* U, double* S, double* vT){

      // wantqo == 0 : as n >= m, only compute S; 
      // wantqo == 1 : as n < m, only compute S; 
      // wantqo == 2 : as n >= m, only compute U, S; 
      // wantqo == 3 : as n < m, only compute U, S; 
      // wantqo == 4 : as n >= m, only compute S, vT; 
      // wantqo == 5 : as n < m, only compute S, vT; 
      // wantqo == 6 : as n >= m, only compute U, S, vT; 
      // wantqo == 7 : as n < m, only compute U, S, vT; 

      // The rules of counting wantqo idx;
      uni10_int wantqo = 0;

      if(!(U == NULL))
        wantqo += 1;

      if(!(vT == NULL))
        wantqo += 2;

      wantqo = (N<M) ? 2 * wantqo + 1 : 2 * wantqo;

      uni10_int min = std::min(M, N);
      uni10_int ldA = N, ldu = N, ldvT = min;
      uni10_int lwork = -1;
      double worktest;
      uni10_int info;
      uni10_int* iwork = (uni10_int*)malloc(8*min*sizeof(uni10_int));

      double *Mij = NULL;
      double *meta_U = U, *meta_vT = vT; 

      char jobz[1];

      //pruni10_intf("wantqo: %d\n", wantqo);

      if(wantqo == 0 || wantqo == 1){
        jobz[0] = 'N';
        Mij = (double*)malloc(M*N*sizeof(double));
        memcpy(Mij, Mij_ori, M*N*sizeof(double));
      }

      else if(wantqo == 2){
        // wantqo == 2 : as n >= m, only compute U, S
        // U --> U && vT --> Mij
        // Where U has been initialized && vT is NULL;
        jobz[0] = 'O';
        Mij = (double*)malloc(M*N*sizeof(double));
        memcpy(Mij, Mij_ori, M*N*sizeof(double));
        meta_vT = NULL;
      }

      else if(wantqo == 3){
        // wantqo == 3 : as n < m, only compute U, S; 
        // vT --> vT && U --> Mij
        // Where U has been initialized && vT is NULL;
        jobz[0] = 'O';
        memcpy(U, Mij_ori, M*N*sizeof(double));
        Mij = U;
        meta_U = NULL;
        meta_vT = (double*)malloc(N*N*sizeof(double));
      }

      else if(wantqo == 4){
        // wantqo == 4 : as n >= m, only compute S, vT; 
        // U --> U && vT --> Mij
        // Where U is NULL && vT has been initialized;
        jobz[0] = 'O';
        memcpy(vT, Mij_ori, M*N*sizeof(double));
        Mij = vT;
        meta_U = (double*)malloc(M*M*sizeof(double));
        meta_vT = NULL;
      }

      else if(wantqo == 5){
        // wantqo == 5 : as n < m, only compute S, vT; 
        // vT --> vT && U --> Mij
        // Where U is NULL && vT has been initialized;
        jobz[0] = 'O';
        Mij = (double*)malloc(M*N*sizeof(double));
        memcpy(Mij, Mij_ori, M*N*sizeof(double));
        meta_U = NULL;
      }

      else if(wantqo == 6){
        // wantqo == 6 : as n >= m, only compute U, S, vT; 
        // U --> U && vT --> Mij
        // Where U && vT have been initialized;
        jobz[0] = 'O';
        memcpy(vT, Mij_ori, M*N*sizeof(double));
        Mij = vT;
        meta_vT = NULL;
      }

      else if(wantqo == 7){
        // wantqo == 7 : as n < m, only compute U, S, vT; 
        // vT --> vT && U --> Mij
        jobz[0] = 'O';
        memcpy(U, Mij_ori, M*N*sizeof(double));
        Mij = U;
        meta_U = NULL;

      }

      dgesdd(jobz, &N, &M, Mij, &ldA, S, meta_vT, &ldu, meta_U, &ldvT, &worktest, &lwork, iwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgesdd': Lapack INFO = ", info);

      lwork = (uni10_int)worktest;
      double *work = (double*)malloc(lwork*sizeof(double));
      dgesdd(jobz, &N, &M, Mij, &ldA, S, meta_vT, &ldu, meta_U, &ldvT, work, &lwork, iwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgesdd': Lapack INFO = ", info);

      if(wantqo == 0 || wantqo == 1 || wantqo == 2 || wantqo == 5)
        free(Mij);

      else if(wantqo == 3)
        free(meta_vT);

      else if(wantqo == 4)
        free(meta_U);

      free(work);
      free(iwork);

    }

    // lapack is builded by fortran which is load by column, so we use 
    // dorgqr -> lq
    // dorglq -> qr
    // dorgrq -> ql 
    // dorgql -> rq
    void matrixQR_cpu(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* R){

      uni10_error_msg(M < N, "%s", "M must be larger than N in matrixQR()");

      memcpy(Q, Mij_ori, N*M*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      uni10_int lda = N;
      uni10_int lwork = -1;
      double worktestdge;
      double worktestdor;
      uni10_int info;
      uni10_int K = N;
      dgelqf(&N, &M, Q, &lda, tau, &worktestdge, &lwork, &info);
      dorglq(&N, &M, &K, Q, &lda, tau, &worktestdor, &lwork, &info);
      lwork = (uni10_int)worktestdge;
      double* workdge = (double*)malloc(lwork*sizeof(double));
      dgelqf(&N, &M, Q, &lda, tau, workdge, &lwork, &info);
      //getR
      uni10_getUpTri(Q, R, M, N);
      //getQ
      lwork = (uni10_int)worktestdor;
      double* workdor = (double*)malloc(lwork*sizeof(double));
      dorglq(&N, &M, &K, Q, &lda, tau, workdor, &lwork, &info);

      //double alpha = 1, beta = 0;
      //dgemm((char*)"N", (char*)"T", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

      free(tau);
      free(workdge);
      free(workdor);
    }

    void matrixRQ_cpu(double* Mij_ori, uni10_int M, uni10_int N, double* R, double* Q){

      uni10_error_msg(N < M, "%s", "N must be larger than M in matrixRQ()");

      memcpy(Q, Mij_ori, M*N*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      uni10_int lda = N;
      uni10_int lwork = -1;
      double worktestdge;
      double worktestdor;
      uni10_int info;
      uni10_int K = M;
      dgeqlf(&N, &M, Q, &lda, tau, &worktestdge, &lwork, &info);
      dorgql(&N, &M, &K, Q, &lda, tau, &worktestdor, &lwork, &info);
      lwork = (uni10_int)worktestdge;
      double* workdge = (double*)malloc(lwork*sizeof(double));
      dgeqlf(&N, &M, Q, &lda, tau, workdge, &lwork, &info);
      //getR
      uni10_getUpTri(Q, R, M, N);
      ///getQ
      lwork = (uni10_int)worktestdor;
      double* workdor = (double*)malloc(lwork*sizeof(double));
      dorgql(&N, &M, &K, Q, &lda, tau, workdor, &lwork, &info);

      //double alpha = 1, beta = 0;
      //dgemm((char*)"T", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, R, &M);

      free(tau);
      free(workdge);
      free(workdor);

    }

    void matrixLQ_cpu(double* Mij_ori, uni10_int M, uni10_int N, double* L, double* Q){

      uni10_error_msg(N < M, "%s","N must be larger than M in matrixLQ()");

      memcpy(Q, Mij_ori, M*N*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      uni10_int lda = N;
      uni10_int lwork = -1;
      double worktestdge;
      double worktestdor;
      uni10_int info;
      uni10_int K = M;
      dgeqrf(&N, &M, Q, &lda, tau, &worktestdge, &lwork, &info);
      dorgqr(&N, &M, &K, Q, &lda, tau, &worktestdor, &lwork, &info);
      lwork = (uni10_int)worktestdge;
      double* workdge = (double*)malloc(lwork*sizeof(double));
      dgeqrf(&N, &M, Q, &lda, tau, workdge, &lwork, &info);
      //getL
      uni10_getDnTri(Q, L, M, N);
      //getQ
      lwork = (uni10_int)worktestdor;
      double* workdor = (double*)malloc(lwork*sizeof(double));
      dorgqr(&N, &M, &K, Q, &lda, tau, workdor, &lwork, &info);

      //double alpha = 1, beta = 0;
      //dgemm((char*)"T", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, L, &M);

      free(tau);
      free(workdge);
      free(workdor);
    }

    void matrixQL_cpu(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* L){

      uni10_error_msg(M < N, "%s", "M must be larger than N in matrixQL()");

      memcpy(Q, Mij_ori, N*M*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      uni10_int lda = N;
      uni10_int lwork = -1;
      double worktestdge;
      double worktestdor;
      uni10_int info;
      uni10_int K = N;
      dgerqf(&N, &M, Q, &lda, tau, &worktestdge, &lwork, &info);
      dorgrq(&N, &M, &K, Q, &lda, tau, &worktestdor, &lwork, &info);
      lwork = (uni10_int)worktestdge;
      double* workdge = (double*)malloc(lwork*sizeof(double));
      dgerqf(&N, &M, Q, &lda, tau, workdge, &lwork, &info);
      //getR
      uni10_getDnTri(Q, L, M, N);
      //getQ
      lwork = (uni10_int)worktestdor;
      double* workdor = (double*)malloc(lwork*sizeof(double));
      dorgrq(&N, &M, &K, Q, &lda, tau, workdor, &lwork, &info);

      //double alpha = 1, beta = 0;
      //dgemm((char*)"N", (char*)"T", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

      free(tau);
      free(workdge);
      free(workdor);
    }

    void matrixQDR_cpu(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* D, double* R){

      uni10_error_msg(M < N, "%s", "M must be larger than N in matrixQDR()");

      memcpy(Q, Mij_ori, N*M*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      uni10_int lda = N;
      uni10_int lwork = -1;
      double worktestdge;
      double worktestdor;
      uni10_int info;
      uni10_int K = N;
      dgelqf(&N, &M, Q, &lda, tau, &worktestdge, &lwork, &info);
      dorglq(&N, &M, &K, Q, &lda, tau, &worktestdor, &lwork, &info);
      lwork = (uni10_int)worktestdge;
      double* workdge = (double*)malloc(lwork*sizeof(double));
      dgelqf(&N, &M, Q, &lda, tau, workdge, &lwork, &info);
      //getR
      uni10_getUpTri(Q, R, M, N);
      uni10_getDiag(R, D, N, N, N);
      for(uni10_int i = 0; i < N; i++)
        for(uni10_int j = 0; j < N-i; j++)
          R[i*N+i+j] /= D[i];
      //getQ
      lwork = (uni10_int)worktestdor;
      double* workdor = (double*)malloc(lwork*sizeof(double));
      dorglq(&N, &M, &K, Q, &lda, tau, workdor, &lwork, &info);

      //double alpha = 1, beta = 0;
      //dgemm((char*)"N", (char*)"T", &N, &N, &M, &alpha, Mij_ori, &N, Mij, &N, &beta, R, &N);

      free(tau);
      free(workdge);
      free(workdor);
    }

    void matrixQDRCPIVOT_cpu(double* Mij_ori, uni10_int M, uni10_int N, double* Q, double* D, double* R){
      uni10_error_msg(M != N, "%s" ,"lack of error");     // D(1 x M): pivoted diagonals
      double* Mij = (double*)malloc(M * N * sizeof(double)); // Q(M x M): orthogonal basis
      setTranspose_cpu(Mij_ori, M, N, Mij); //column major + square matrix  // TP(M x M)
      uni10_int max = M > N ? M : N;
      uni10_int min = M < N ? M : N;
      uni10_int lda = max;
      uni10_int lwork = -1;
      double worktest;
      uni10_int* jpvt = (uni10_int*)malloc(N * sizeof(uni10_int));	//column vectors
      memset(jpvt, 0, N * sizeof(uni10_int));
      double* tau = (double*)malloc(min * sizeof(double));
      uni10_int info;
      dgeqp3(&M, &N, Mij, &lda, jpvt, tau, &worktest, &lwork, &info);
      lwork = (uni10_int)worktest;
      double *work = (double*) malloc(lwork * sizeof(double));
      dgeqp3(&M, &N, Mij, &lda, jpvt, tau, work, &lwork, &info);
      uni10_error_msg(info != 0, "%s" ,"lack of error");     // D(1 x M): pivoted diagonals
      for(uni10_int i = 0; i < M; i++)
        D[i] = Mij[i * N + i];  // D Here!!!
      double* T = (double*)calloc(M * N, sizeof(double));
      for(uni10_int i = 0; i < M; i++)
        for(uni10_int j = i; j < N; j++)
          if(i == j)
            T[i * N + j] = 1;
          else
            T[i * N + j] = Mij[j * N + i] / D[i];
      for(uni10_int i = 0; i < M; i++)
        for(uni10_int j = 0; j < N; j++)
          R[i * N + (jpvt[j]-1)] = T[i * N + j];	// TP Here!!!
      dorgqr(&M, &N, &N, Mij, &lda, tau, work, &lwork, &info);
      uni10_error_msg(info != 0, "%s" ,"lack of error");     // D(1 x M): pivoted diagonals
      setTranspose_cpu(Mij, M, N, Q);	// Q Here!!!
      free(Mij);
      free(work);
      free(jpvt);
      free(tau);
      free(T);

    }

    void matrixLDQ_cpu(double* Mij_ori, uni10_int M, uni10_int N, double* L, double* D, double* Q){

      uni10_error_msg(N < M, "%s","N must be larger than M in matrixLDQ()");

      memcpy(Q, Mij_ori, M*N*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      uni10_int lda = N;
      uni10_int lwork = -1;
      double worktestdge;
      double worktestdor;
      uni10_int info;
      uni10_int K = M;
      dgeqrf(&N, &M, Q, &lda, tau, &worktestdge, &lwork, &info);
      dorgqr(&N, &M, &K, Q, &lda, tau, &worktestdor, &lwork, &info);
      lwork = (uni10_int)worktestdge;
      double* workdge = (double*)malloc(lwork*sizeof(double));
      dgeqrf(&N, &M, Q, &lda, tau, workdge, &lwork, &info);
      //get D and L
      uni10_getDnTri(Q, L, M, N);
      uni10_getDiag(L, D, M, M, M);
      for(uni10_int i = 0; i < M; i++)
        for(uni10_int j = 0; j < M-i; j++)
          L[(i+j)*M+i] /= D[i];
      //getQ
      lwork = (uni10_int)worktestdor;
      double* workdor = (double*)malloc(lwork*sizeof(double));
      dorgqr(&N, &M, &K, Q, &lda, tau, workdor, &lwork, &info);
      //getR
      //double alpha = 1, beta = 0;
      //dgemm((char*)"T", (char*)"N", &M, &M, &N, &alpha, Mij, &N, Mij_ori, &N, &beta, L, &M);

      free(tau);
      free(workdge);
      free(workdor);
    }

    void matrixInv_cpu(double* A, uni10_int N){

      uni10_int *ipiv = (uni10_int*)malloc((N+1)*sizeof(uni10_int));
      uni10_int info;
      dgetrf(&N, &N, A, &N, ipiv, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgetrf': Lapack INFO = ", info);

      uni10_int lwork = -1;
      double worktest = 0.;
      dgetri(&N, A, &N, ipiv, &worktest, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgetri': Lapack INFO = ", info);

      lwork = (uni10_int)worktest;
      double *work = (double*)malloc(lwork * sizeof(double));
      dgetri(&N, A, &N, ipiv, work, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dgetri': Lapack INFO = ", info);

      free(ipiv);
      free(work);
    }

    double matrixDet_cpu(double* A, uni10_int N){

      uni10_int *ipiv = (uni10_int*)malloc((N+1)*sizeof(uni10_int));
      uni10_int lwork = 64 * N;
      double *work = (double*)malloc(lwork * sizeof(double));
      uni10_int info;
      dgetrf(&N,&N,A,&N,ipiv,&info);
      uni10_error_msg( info != 0, "%s %d", "Error in Lapack function 'dgetrf': Lapack INFO = ", info );
      double det = 1;
      uni10_int neg = 0;
      for (uni10_int i = 0; i < N; i++) {
        det *= A[i * N + i];
        if (ipiv[i] != (i+1)) neg = !neg;
      }
      free(ipiv);
      free(work);
      return neg?-det:det;

    }

    void setIdentity_cpu(double* elem, uni10_uint64 M, uni10_uint64 N){
      uni10_uint64 min;
      if(M < N)
        min = M;
      else  
        min = N;
      memset(elem, 0, M * N * sizeof(double));
      for(uni10_uint64 i = 0; i < min; i++)
        elem[i * N + i] = 1;
    }

    void trimatrixEigH_cpu(double* D, double* E, uni10_int N, 
        double* z, uni10_int LDZ){

      uni10_int info;
      double* work;

      if(z==NULL){
        work = NULL;
        LDZ = N;
        dstev((char*)"N", &N, D, E, z, &LDZ, work, &info);
      }else{
        work         = (double*)malloc(4*N*sizeof(double));
        dstev((char*)"V", &N, D, E, z, &LDZ, work, &info);
        free(work);
      }

    }

    void eigSyDecompose_cpu(double* Kij, uni10_int N, double* Eig, double* EigVec){

      memcpy(EigVec, Kij, N * N * sizeof(double));
      uni10_int ldA = N;
      uni10_int lwork = -1;
      double worktest;
      uni10_int info;
      dsyev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, &worktest, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dsyev': Lapack INFO = ", info);

      lwork = (uni10_int)worktest;
      double* work= (double*)malloc(sizeof(double)*lwork);
      dsyev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, work, &lwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'dsyev': Lapack INFO = ", info);

      free(work);
    }

    void matrix_diag_dense_add_cpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b){

      memcpy(b, a, m*n*sizeof(double));
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }

    void matrix_diag_dense_sub_cpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* b){

      uni10_uint64 elemNum = m*n;
      memcpy(b, a, elemNum*sizeof(uni10_double64));
      vectorScal_cpu(-1., b, elemNum);

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }

    void matrix_diag_dense_mul_cpu(const double* D, const double* a, uni10_uint64 m, uni10_uint64 n, double* v){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        v[i] = D[i] * a[i*n+i];

    }

    void matrix_dense_diag_add_cpu(double* a, const double* D, uni10_uint64 m, uni10_uint64 n){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        a[i*n+i] += D[i];

    }
                                                                                         
    void matrix_dense_diag_sub_cpu(double* a, const double* D, uni10_uint64 m, uni10_uint64 n){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        a[i*n+i] -= D[i];

    }

    void matrix_diag_dense_mul_cpu(double* D, const double* a, uni10_uint64 m, uni10_uint64 n){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        D[i] *= a[i*n+i];

    }
   
    void matrix_dense_diag_add_cpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b){

      memcpy(b, a, m*n*sizeof(double));
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }

    void matrix_dense_diag_sub_cpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* b){

      memcpy(b, a, m*n*sizeof(double));
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }

    void matrix_dense_diag_mul_cpu(const double* a, const double* D, uni10_uint64 m, uni10_uint64 n, double* v){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        v[i] = a[i*n+i] * D[i];

    }

  } /* namespace uni10_linalg */

} /* namespace uni10 */

