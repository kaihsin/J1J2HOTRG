#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_linalg_lapack_cpu_dz.h"


namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(double a, double* X, uni10_int incx, std::complex<double>* Y, uni10_int incy, uni10_uint64 N){   // Y = aX + Y

      double* cr = (double*)Y;

      int64_t left      = N;
      uni10_uint64 offset_r = 0;
      uni10_uint64 offset_c = 0;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        daxpy(&chunk, &a, X + offset_r, &incx, cr + offset_c, &incy);
        offset_r += chunk;
        offset_c += 2 * (uni10_uint64)chunk;
        left -= INT_MAX;
      }

    }

    void vectorAdd(std::complex<double>* Y, double* X, uni10_uint64 N){

      double* cr = (double*)Y;
      uni10_double64 a    = 1.0;

      int64_t left        = N;
      uni10_int incx      = 1;
      uni10_int incy      = 2;
      uni10_uint64 offset_r = 0;
      uni10_uint64 offset_c = 0;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        daxpy(&chunk, &a, X + offset_r, &incx, cr + offset_c, &incy);
        offset_r += chunk;
        offset_c += 2 * (uni10_uint64)chunk;
        left -= INT_MAX;
      }

    }

    void vectorSub(std::complex<double>* Y, double* X, uni10_uint64 N){

      double* cr = (double*)Y;
      uni10_double64 a    = -1.0;

      int64_t left        = N;
      uni10_int incx      = 1;
      uni10_int incy      = 2;
      uni10_uint64 offset_r = 0;
      uni10_uint64 offset_c = 0;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        daxpy(&chunk, &a, X + offset_r, &incx, cr + offset_c, &incy);
        offset_r += chunk;
        offset_c += 2 * (uni10_uint64)chunk;
        left -= INT_MAX;
      }

    }

    void vectorMul(std::complex<double>* Y, double* X, uni10_uint64 N){

      for(uni10_uint64 i = 0; i < N; i++)
        Y[i] *= X[i];

    }

    void vectorScal(double a, std::complex<double>* X, uni10_uint64 N){

      int64_t left = N;
      uni10_int inc = 1;
      uni10_uint64 offset = 0;
      uni10_int chunk;

      while(left > 0){
        if(left > INT_MAX)
          chunk = INT_MAX;
        else
          chunk = left;
        zdscal(&chunk, &a, X + offset, &inc);
        offset += chunk;
        left -= INT_MAX;
      }

    }

    void vectorExp(double a, std::complex<double>* X, uni10_uint64 N){
      for(uni10_uint64 i = 0; i < N; i++)
        X[i] = std::exp(a * X[i]);
    }

    void matrixDot(double* A, std::complex<double>* B, uni10_int M, uni10_int N, uni10_int K, std::complex<double>* C){

      uni10_int size_A = M * K;
      std::complex<double>* CA = (std::complex<double>*)malloc(size_A*sizeof(std::complex<double>));
      uni10_elem_cast(CA, A, size_A);
      std::complex<double> alpha = 1.0, beta = 0.0;
      zgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, B, &N, CA, &K, &beta, C, &N);
      free(CA);
    }

    void matrixDot(std::complex<double>* A, double* B, uni10_int M, uni10_int N, uni10_int K, std::complex<double>* C){

      uni10_int size_B = K * N;
      std::complex<double>* CB = (std::complex<double>*)malloc(size_B*sizeof(std::complex<double>));
      std::complex<double> alpha = 1.0, beta = 0.0;
      uni10_elem_cast(CB, B, size_B);
      zgemm((char*)"N", (char*)"N", &N, &M, &K, &alpha, CB, &N, A, &K, &beta, C, &N);
      free(CB);
    }

    void diagRowMul(std::complex<double>* mat, double* diag, uni10_uint64 M, uni10_uint64 N){
      for(uni10_uint64 i = 0; i < M; i++)
        vectorScal(diag[i], &(mat[i * N]), N);
    }

    void diagColMul(std::complex<double>* mat, double* diag, uni10_uint64 M, uni10_uint64 N){
      for(uni10_uint64 i = 0; i < M; i++){
        uni10_uint64 ridx = i * N;
        for(uni10_uint64 j = 0; j < N; j++)
          mat[ridx + j] *= diag[j];
      }
    }

    void eigDecompose(double* Kij_ori, uni10_int N, std::complex<double>* Eig, std::complex<double>* EigVec){
      std::complex<double> *Kij = (std::complex<double>*) malloc(N * N * sizeof(std::complex<double>));
      uni10_elem_cast(Kij, Kij_ori, N * N);
      eigDecompose(Kij, N, Eig, EigVec);
      free(Kij);
    }

    void eigSyDecompose(std::complex<double>* Kij, uni10_int N, double* Eig, std::complex<double>* EigVec){
      //eigDecompose(Kij, N, Eig, EigVec, ongpu);
      memcpy(EigVec, Kij, N * N * sizeof(std::complex<double>));
      uni10_int ldA = N;
      uni10_int lwork = -1;
      std::complex<double> worktest;
      double* rwork = (double*) malloc((3*N+1) * sizeof(double));
      uni10_int info;
      zheev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, &worktest, &lwork, rwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zheev': Lapack INFO = ", info);

      lwork = (uni10_int)worktest.real();
      std::complex<double>* work= (std::complex<double>*)malloc(sizeof(std::complex<double>)*lwork);
      zheev((char*)"V", (char*)"U", &N, EigVec, &ldA, Eig, work, &lwork, rwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zheev': Lapack INFO = ", info);

      free(work);
      free(rwork);
    }

    void matrixSVD(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, double *S, std::complex<double>* vT){

      char jobu[1], jobv[1]; 
      jobu[0] = ( U  == NULL ) ? 'N' : 'S';
      jobv[0] = ( vT == NULL ) ? 'N' : 'S';

      std::complex<double>* Mij = (std::complex<double>*)malloc(M * N * sizeof(std::complex<double>));
      memcpy(Mij, Mij_ori, M * N * sizeof(std::complex<double>));
      uni10_int min = std::min(M, N);
      uni10_int ldA = N, ldu = N, ldvT = min;
      uni10_int lwork = -1;
      std::complex<double> worktest;
      uni10_int info;
      double *rwork = (double*) malloc(std::max( (uni10_int)1, 5 * min) * sizeof(double));
      zgesvd(jobv, jobu, &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, &worktest, &lwork, rwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zgesvd': Lapack INFO = ", info);

      lwork = (uni10_int)(worktest.real());
      std::complex<double> *work = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgesvd(jobv, jobu, &N, &M, Mij, &ldA, S, vT, &ldu, U, &ldvT, work, &lwork, rwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zgesvd': Lapack INFO = ", info);

      free(rwork);
      free(work);
      free(Mij);
    }

    void matrixSDD(std::complex<double>* Mij_ori, uni10_int M, uni10_int N, std::complex<double>* U, double *S, std::complex<double>* vT){

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
      std::complex<double> worktest;
      uni10_int info;
      uni10_int* iwork = (uni10_int*)malloc(8*min*sizeof(uni10_int));

      std::complex<double> *Mij = NULL;
      std::complex<double> *meta_U = U, *meta_vT = vT; 
      double* rwork = (double*)malloc((5*min*min + 7 * min) * sizeof(double));

      char jobz[1];

      //pruni10_intf("wantqo: %d\n", wantqo);

      if(wantqo == 0 || wantqo == 1){
        jobz[0] = 'N';
        Mij = (std::complex<double>*)malloc(M*N*sizeof(std::complex<double>));
        memcpy(Mij, Mij_ori, M*N*sizeof(std::complex<double>));
      }

      else if(wantqo == 2){
        // wantqo == 2 : as n >= m, only compute U, S
        // U --> U && vT --> Mij
        // Where U has been initialized && vT is NULL;
        jobz[0] = 'O';
        Mij = (std::complex<double>*)malloc(M*N*sizeof(std::complex<double>));
        memcpy(Mij, Mij_ori, M*N*sizeof(std::complex<double>));
        meta_vT = NULL;
      }

      else if(wantqo == 3){
        // wantqo == 3 : as n < m, only compute U, S; 
        // vT --> vT && U --> Mij
        // Where U has been initialized && vT is NULL;
        jobz[0] = 'O';
        memcpy(U, Mij_ori, M*N*sizeof(std::complex<double>));
        Mij = U;
        meta_U = NULL;
        meta_vT = (std::complex<double>*)malloc(N*N*sizeof(std::complex<double>));
      }

      else if(wantqo == 4){
        // wantqo == 4 : as n >= m, only compute S, vT; 
        // U --> U && vT --> Mij
        // Where U is NULL && vT has been initialized;
        jobz[0] = 'O';
        memcpy(vT, Mij_ori, M*N*sizeof(std::complex<double>));
        Mij = vT;
        meta_U = (std::complex<double>*)malloc(M*M*sizeof(std::complex<double>));
        meta_vT = NULL;
      }

      else if(wantqo == 5){
        // wantqo == 5 : as n < m, only compute S, vT; 
        // vT --> vT && U --> Mij
        // Where U is NULL && vT has been initialized;
        jobz[0] = 'O';
        Mij = (std::complex<double>*)malloc(M*N*sizeof(std::complex<double>));
        memcpy(Mij, Mij_ori, M*N*sizeof(std::complex<double>));
        meta_U = NULL;
      }

      else if(wantqo == 6){
        // wantqo == 6 : as n >= m, only compute U, S, vT; 
        // U --> U && vT --> Mij
        // Where U && vT have been initialized;
        jobz[0] = 'O';
        memcpy(vT, Mij_ori, M*N*sizeof(std::complex<double>));
        Mij = vT;
        meta_vT = NULL;
      }

      else if(wantqo == 7){
        // wantqo == 7 : as n < m, only compute U, S, vT; 
        // vT --> vT && U --> Mij
        jobz[0] = 'O';
        memcpy(U, Mij_ori, M*N*sizeof(std::complex<double>));
        Mij = U;
        meta_U = NULL;

      }

      zgesdd(jobz, &N, &M, Mij, &ldA, S, meta_vT, &ldu, meta_U, &ldvT, &worktest, &lwork, rwork, iwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zgesdd': Lapack INFO = ", info);

      lwork = (uni10_int)(worktest.real());
      std::complex<double> *work = (std::complex<double>*)malloc(lwork*sizeof(std::complex<double>));
      zgesdd(jobz, &N, &M, Mij, &ldA, S, meta_vT, &ldu, meta_U, &ldvT, work, &lwork, rwork, iwork, &info);

      uni10_error_msg(info != 0, "%s %d", "Error in Lapack function 'zgesdd': Lapack INFO = ", info);

      if(wantqo == 0 || wantqo == 1 || wantqo == 2 || wantqo == 5)
        free(Mij);

      else if(wantqo == 3)
        free(meta_vT);

      else if(wantqo == 4)
        free(meta_U);

      free(work);
      free(rwork);
      free(iwork);

    }

    // r operator() z

    void matrix_diag_dense_add(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      memcpy(b, a, m*n*sizeof(uni10_complex128));
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }

    void matrix_diag_dense_sub(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_uint64 elemNum = m*n;
      memcpy(b, a, elemNum*sizeof(uni10_complex128));
      vectorScal(-1., b,elemNum);
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }

    void matrix_diag_dense_mul(const double* D, const std::complex<double>* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* v){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        v[i] = D[i] - a[i*n+i];

    }
   
    void matrix_dense_diag_add(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_elem_cast(b, a, m*n);
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];


    }
                                                                                                                                      
    void matrix_dense_diag_sub(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_elem_cast(b, a, m*n);
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] -= D[i];

    }
                                                                                                                                      
    void matrix_dense_diag_mul(const double* a, const std::complex<double>* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* v){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        v[i] = D[i] * a[i*n+i];

    }

    // z operator() r

    void matrix_diag_dense_add(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_elem_cast(b, a, m*n);
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }

    void matrix_diag_dense_sub(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_uint64 elemNum = m*n;
      uni10_elem_cast(b, a, elemNum);
      vectorScal(-1., b, elemNum);

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }

    void matrix_diag_dense_mul(const std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n, std::complex<double>* v){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        v[i] = D[i] * a[i*n+i];

    }

    void matrix_dense_diag_add(std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        a[i*n+i] += D[i];

    }
                                                                                                       
    void matrix_dense_diag_sub(std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        a[i*n+i] -= D[i];

    }
                                                                                                       
    void matrix_diag_dense_mul(std::complex<double>* D, const double* a, uni10_uint64 m, uni10_uint64 n){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        D[i] *= a[i*n+i];

    }
   
    void matrix_dense_diag_add(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_elem_copy(b, a, m*n*sizeof(uni10_complex128));
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] += D[i];

    }
                                                                                                                                       
    void matrix_dense_diag_sub(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* b){

      uni10_elem_copy(b, a, m*n*sizeof(uni10_complex128));
      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        b[i*n+i] -= D[i];

    }
                                                                                                                                       
    void matrix_dense_diag_mul(const std::complex<double>* a, const double* D, uni10_uint64 m, uni10_uint64 n, std::complex<double>* v){

      uni10_uint64 min = std::min(m, n);
      uni10_uint64 i;
      for( i = 0; i < min; i++)
        v[i] = D[i] * a[i*n+i];

    }

  };/* namespace uni10_linalg */

};/* namespace uni10 */

