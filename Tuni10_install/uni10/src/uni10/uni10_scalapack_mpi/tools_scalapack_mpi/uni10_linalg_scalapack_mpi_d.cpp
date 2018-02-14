#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_tools_scalapack_mpi.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_linalg_scalapack_mpi_d.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_blacs_wrapper_mpi.h"

#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_scalapack_wrapper_mpi.h"

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(uni10_double64 a, uni10_double64* X, int incx, uni10_double64* Y, int incy, size_t N){   // Y = aX + Y

      uni10_uint64 left = N;
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

      MPI_Barrier(MPI_COMM_WORLD);

    }

    void vectorAdd(uni10_double64* Y, uni10_double64* X, uni10_uint64 N){   // Y = Y + X

      uni10_double64 a    = 1.0;
      uni10_int inc       = 1;
      uni10_int left      = N;
      uni10_uint64 offset = 0;
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

      MPI_Barrier(MPI_COMM_WORLD);

    }

    void vectorSub(uni10_double64* Y, uni10_double64* X, size_t N){   // Y = Y + X

      uni10_double64 a    = -1.0;
      uni10_int inc       = 1;
      uni10_int left      = N;
      uni10_uint64 offset = 0;
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

      MPI_Barrier(MPI_COMM_WORLD);

    }

    void vectorMul(uni10_double64* Y, uni10_double64* X, uni10_uint64 N){ // Y = Y * X, element-wise multiplication;
      
      for(uni10_uint64 i = 0; i < N; i++)
        Y[i] *= X[i];

    }

    void vectorScal(uni10_double64 a, uni10_double64* X, uni10_uint64 N){

      uni10_int inc = 1;
      uni10_int left = N;
      uni10_uint64 offset = 0;
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

    void vectorExp(uni10_double64 a, uni10_double64* X, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    uni10_double64 vectorSum(const uni10_double64* X, uni10_uint64 mpxnq, uni10_int inc){

      uni10_int rank = env_variables.get_info().rank_mpi;
      uni10_double64 sum = 0;
      uni10_double64 subsum = 0;
      uni10_uint64 idx = 0;
      for(uni10_uint64 i = 0; i < mpxnq; i++){
        subsum += X[idx];
        idx += inc;
      }

      MPI_Allreduce(&subsum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      return sum;
      
    }

    uni10_double64 vectorNorm(const uni10_double64* A, const uni10_int* descA, uni10_uint64 M, uni10_uint64 N){
      
      uni10_int ione = 1;
      return pdlange("F", &N, &M, A, &ione, &ione, descA, NULL);
    }

    void matrixDot(const uni10_double64* A, const uni10_int* descA, const uni10_double64* B, const uni10_int* descB, 
        const uni10_int M, const uni10_int N, const uni10_int K, uni10_double64* C, uni10_int* descC){

      uni10_double64 alpha = 1.0, beta = 0.;
      uni10_int ione      = 1;
      pdgemm("N", "N", &N, &M, &K, &alpha, B, &ione, &ione, descB, A, &ione, &ione, descA, 
          &beta, C, &ione, &ione, descC);

    }

    void diagRowMul(uni10_double64* mat, uni10_double64* diag, uni10_int mp, uni10_int nq, uni10_int ch){
      
      for(uni10_int i = 0; i < nq; i++){
        uni10_int ridx = i * mp;
        for(uni10_int j = 0; j < mp; j++)
          mat[ridx + j] *= diag[ch+i];
      }

    }

    void diagColMul(uni10_double64 *mat, uni10_double64* diag, uni10_int mp, uni10_int nq, uni10_int rh){

      for(uni10_int i = 0; i < nq; i++){
        uni10_uint64 ridx = i * mp;
        for(uni10_int j = 0; j < mp; j++)
          mat[ridx + j] *= diag[rh+j];
      }

    }

    void setTranspose(const uni10_double64* A, const uni10_int* descA, uni10_int M, uni10_int N, uni10_double64* AT, uni10_int* descAT){

      uni10_double64 alpha = 1.0;
      uni10_double64 beta  = 0.;
      uni10_int ione = 1;
      uni10_const_bool balance = env_variables.get_info().nprow == env_variables.get_info().npcol;

      if(balance){

        pdgeadd_((char*)"T", &N, &M, &alpha, A, &ione, &ione, descA, &beta, AT, &ione, &ione, descAT);

      }else{

        uni10_error_msg(!balance, "%s", "Still dosen't support if nprow != npcol. [Developping]");

      }

    }

    void setDagger(const uni10_double64* A, const uni10_int* descA, uni10_int M, uni10_int N, uni10_double64* AT, uni10_int* descAT){

      setTranspose(A, descA, M, N, AT, descAT);

    }

    void matrixSVD(const uni10_double64* Mij_ori, const uni10_int mp, const uni10_int nq, const uni10_int* descMij_ori, 
        uni10_int M, uni10_int N, uni10_double64* U, uni10_int* descU, uni10_double64* S, uni10_double64* vT, uni10_int* descvT){

      char jobu[1], jobv[1]; 
      jobu[0] = ( descU  == NULL ) ? 'N' : 'V';
      jobv[0] = ( descvT == NULL ) ? 'N' : 'V';

      uni10_int descMij[9];
      memcpy(descMij, descMij_ori, 9*sizeof(uni10_int));

      uni10_int memsize = mp * nq * sizeof(uni10_double64);
      uni10_double64* Mij = (uni10_double64*)malloc( memsize );
      memcpy(Mij, Mij_ori, memsize);


      uni10_int mb = env_variables.get_info().blockgrid;
      uni10_int nb = env_variables.get_info().blockgrid;

      printDesc(descMij);
      printDesc(descU);
      printDesc(descvT);
      exit(0);

      //uni10_int MYROW = env_variables.get_info().myrow;
      //uni10_int MYCOL = env_variables.get_info().mycol;
      //uni10_int NPROW = env_variables.get_info().nprow;
      //uni10_int NPCOL = env_variables.get_info().npcol;

      //uni10_int mp1 = numroc(&N, &mb, &MYROW, &descMij[1], &NPROW);
      //uni10_int nq1 = numroc(&M, &nb, &MYCOL, &descMij[8], &NPCOL);

      //printf("rank: %d, mp: %d, nq: %d, mp1: %d, nq1: %d\n", env_variables.get_info().rank_mpi, mp, nq, mp1, nq1);
      //exit(0);

      uni10_int ione = 1;
      uni10_int info;

      uni10_int min = std::min(M, N);
      uni10_int sizeb = std::max(M, N);

      uni10_int mp0, nq0; 
      if(env_variables.get_info().rank_mpi == 0){
        mp0  = mp;
        nq0  = nq;
      }

      MPI_Bcast(&mp0, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&nq0, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      uni10_int watobd =  nb * ( mp0 + nq0 + 1 ) + nq0;
      
      uni10_int wantu  = ( descU  == NULL ) ? 0 : 1;
      uni10_int wantvt = ( descvT == NULL ) ? 0 : 1;

      uni10_int wbdtosvd = min * (wantu * M + wantvt * N);

      uni10_int bdsqr      =  std::max( 1, 2 * min + ( 2 * min - 4) * std::max(wantu, wantvt));
      uni10_int ormbrqln   =  std::max( (nb*(nb-1))/2, (nq + mp) * nb ) + nb*nb;
      uni10_int ormbrprt   =  std::max( (mb*(mb-1))/2, (mp + nq) * mb ) + mb*mb;

      wbdtosvd += std::max(bdsqr, std::max(wantvt * ormbrqln, wantu * ormbrprt));

      uni10_int lwork = 2 + 6*sizeb + std::max(watobd, wbdtosvd) + 1;

      double *work = (double*)malloc(lwork*sizeof(double));
      MPI_Barrier(MPI_COMM_WORLD);

      pdgesvd(jobv, jobu, &N, &M, Mij, &ione, &ione, descMij, S, vT, &ione, &ione, descvT, U, &ione, &ione, descU, work, &lwork, &info);

      free(Mij);
      free(work);

    }

    void matrixSDD(uni10_double64* Mij_ori, int M, int N, uni10_double64* U, uni10_double64* S, uni10_double64* vT){

      uni10_error_msg(true, "%s", "Don't support in scalapack!!!");

    }

    void matrixQR(const uni10_double64* Mij_ori, const uni10_int mp, const uni10_int nq, const uni10_int* descMij_ori, 
        uni10_int M, uni10_int N, uni10_double64* Q, uni10_int* descQ, uni10_double64* R, uni10_int* descR){

      uni10_error_msg(M < N, "%s", "M must be larger than N in matrixQR()");

      uni10_int ione = 1;
      uni10_int descMij[9];
      memcpy(descMij, descMij_ori, 9*sizeof(uni10_int));

      double* Mij = (double*)malloc(mp*nq*sizeof(double));
      memcpy(Mij, Mij_ori, mp*nq*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      uni10_int info;
      uni10_int mb_a = descMij[4];
      uni10_int nb_a = descMij[5];
      uni10_int rsrc_a = descMij[6];
      uni10_int csrc_a = descMij[7];
      uni10_int MYROW  = env_variables.get_info().myrow;
      uni10_int MYCOL  = env_variables.get_info().mycol;
      uni10_int NPROW  = env_variables.get_info().nprow;
      uni10_int NPCOL  = env_variables.get_info().npcol;
     
      uni10_int iarow = indxg2p(&ione, &mb_a, &MYROW, &rsrc_a, &NPROW);
      uni10_int iacol = indxg2p(&ione, &nb_a, &MYCOL, &csrc_a, &NPCOL);
      uni10_int mp0   = numroc(&N, &mb_a, &MYROW, &iarow, &NPROW);
      uni10_int nq0   = numroc(&M, &nb_a, &MYCOL, &iacol, &NPCOL);
      uni10_int lwork = mb_a * (mp0 + nq0 + mb_a);

      double* work = (double*)malloc(lwork*sizeof(double));

      pdgelqf(&N, &M, Mij, &ione, &ione, descMij, tau, work, &lwork, &info);

      //getR
      pdlacpy("L", &N, &N, Mij, &ione, &ione, descMij, R, &ione, &ione, descR);

      uni10_int K = N;
      pdorglq(&N, &M, &K, Mij, &ione, &ione, descMij, tau, work, &lwork, &info);

      //getQ
      memcpy(Q, Mij, mp*nq*sizeof(double));

      free(Mij);
      free(tau);
      free(work);

    }

    void matrixRQ(const uni10_double64* Mij_ori, const uni10_int mp, const uni10_int nq, const uni10_int* descMij_ori,
        uni10_int M, uni10_int N, uni10_double64* R, uni10_int* descR, uni10_double64* Q, uni10_int* descQ){

      uni10_error_msg(N < M, "%s","N must be larger than M in matrixLQ()");

      uni10_int ione = 1;
      uni10_int descMij[9];
      memcpy(descMij, descMij_ori, 9*sizeof(uni10_int));

      double* Mij = (double*)malloc(mp*nq*sizeof(double));
      memcpy(Mij, Mij_ori, mp*nq*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      uni10_int info;
      uni10_int mb_a = descMij[4];
      uni10_int nb_a = descMij[5];
      uni10_int rsrc_a = descMij[6];
      uni10_int csrc_a = descMij[7];
      uni10_int MYROW  = env_variables.get_info().myrow;
      uni10_int MYCOL  = env_variables.get_info().mycol;
      uni10_int NPROW  = env_variables.get_info().nprow;
      uni10_int NPCOL  = env_variables.get_info().npcol;
     
      uni10_int iarow = indxg2p(&ione, &mb_a, &MYROW, &rsrc_a, &NPROW);
      uni10_int iacol = indxg2p(&ione, &nb_a, &MYCOL, &csrc_a, &NPCOL);
      uni10_int mp0   = numroc(&N, &mb_a, &MYROW, &iarow, &NPROW);
      uni10_int nq0   = numroc(&M, &nb_a, &MYCOL, &iacol, &NPCOL);
      uni10_int lwork = mb_a * (mp0 + nq0 + mb_a);

      double* work = (double*)malloc(lwork*sizeof(double));

      pdgeqlf(&N, &M, Mij, &ione, &ione, descMij, tau, work, &lwork, &info);

      //getR
      pdlacpy("L", &M, &M, Mij, &ione, &ione, descMij, R, &ione, &ione, descR);

      uni10_int K = M;
      pdorgql(&N, &M, &K, Mij, &ione, &ione, descMij, tau, work, &lwork, &info);

      //getQ
      memcpy(Q, Mij, mp*nq*sizeof(double));

      free(Mij);
      free(tau);
      free(work);

    }

    void matrixLQ(const uni10_double64* Mij_ori, const uni10_int mp, const uni10_int nq, const uni10_int* descMij_ori, 
        uni10_int M, uni10_int N, uni10_double64* L, uni10_int* descL, uni10_double64* Q, uni10_int* descQ){

      uni10_error_msg(N < M, "%s","N must be larger than M in matrixLQ()");

      uni10_int ione = 1;
      uni10_int descMij[9];
      memcpy(descMij, descMij_ori, 9*sizeof(uni10_int));

      double* Mij = (double*)malloc(mp*nq*sizeof(double));
      memcpy(Mij, Mij_ori, mp*nq*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      uni10_int info;
      uni10_int mb_a = descMij[4];
      uni10_int nb_a = descMij[5];
      uni10_int rsrc_a = descMij[6];
      uni10_int csrc_a = descMij[7];
      uni10_int MYROW  = env_variables.get_info().myrow;
      uni10_int MYCOL  = env_variables.get_info().mycol;
      uni10_int NPROW  = env_variables.get_info().nprow;
      uni10_int NPCOL  = env_variables.get_info().npcol;
     
      uni10_int iarow = indxg2p(&ione, &mb_a, &MYROW, &rsrc_a, &NPROW);
      uni10_int iacol = indxg2p(&ione, &nb_a, &MYCOL, &csrc_a, &NPCOL);
      uni10_int mp0   = numroc(&N, &mb_a, &MYROW, &iarow, &NPROW);
      uni10_int nq0   = numroc(&M, &nb_a, &MYCOL, &iacol, &NPCOL);
      uni10_int lwork = mb_a * (mp0 + nq0 + mb_a);

      double* work = (double*)malloc(lwork*sizeof(double));

      pdgeqrf(&N, &M, Mij, &ione, &ione, descMij, tau, work, &lwork, &info);

      //getR
      pdlacpy("U", &M, &M, Mij, &ione, &ione, descMij, L, &ione, &ione, descL);

      uni10_int K = M;
      pdorgqr(&N, &M, &K, Mij, &ione, &ione, descMij, tau, work, &lwork, &info);

      //getQ
      memcpy(Q, Mij, mp*nq*sizeof(double));

      free(Mij);
      free(tau);
      free(work);

    }

    void matrixQL(const uni10_double64* Mij_ori, const uni10_int mp, const uni10_int nq, const uni10_int* descMij_ori,
        uni10_int M, uni10_int N, uni10_double64* Q, uni10_int* descQ, uni10_double64* L, uni10_int* descL){

      uni10_error_msg(M < N, "%s", "M must be larger than N in matrixQR()");

      uni10_int ione = 1;
      uni10_int descMij[9];
      memcpy(descMij, descMij_ori, 9*sizeof(uni10_int));

      double* Mij = (double*)malloc(mp*nq*sizeof(double));
      memcpy(Mij, Mij_ori, mp*nq*sizeof(double));
      double* tau = (double*)malloc(M*sizeof(double));
      uni10_int info;
      uni10_int mb_a = descMij[4];
      uni10_int nb_a = descMij[5];
      uni10_int rsrc_a = descMij[6];
      uni10_int csrc_a = descMij[7];
      uni10_int MYROW  = env_variables.get_info().myrow;
      uni10_int MYCOL  = env_variables.get_info().mycol;
      uni10_int NPROW  = env_variables.get_info().nprow;
      uni10_int NPCOL  = env_variables.get_info().npcol;
     
      uni10_int iarow = indxg2p(&ione, &mb_a, &MYROW, &rsrc_a, &NPROW);
      uni10_int iacol = indxg2p(&ione, &nb_a, &MYCOL, &csrc_a, &NPCOL);
      uni10_int mp0   = numroc(&N, &mb_a, &MYROW, &iarow, &NPROW);
      uni10_int nq0   = numroc(&M, &nb_a, &MYCOL, &iacol, &NPCOL);
      uni10_int lwork = mb_a * (mp0 + nq0 + mb_a);

      double* work = (double*)malloc(lwork*sizeof(double));

      pdgerqf(&N, &M, Mij, &ione, &ione, descMij, tau, work, &lwork, &info);

      //getR
      pdlacpy("U", &N, &N, Mij, &ione, &ione, descMij, L, &ione, &ione, descL);

      uni10_int K = N;
      pdorgrq(&N, &M, &K, Mij, &ione, &ione, descMij, tau, work, &lwork, &info);

      //getQ
      memcpy(Q, Mij, mp*nq*sizeof(double));

      free(Mij);
      free(tau);
      free(work);

    }

    void matrixQDR(uni10_double64* Mij_ori, int M, int N, uni10_double64* Q, uni10_double64* D, uni10_double64* R){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixQDRCPIVOT(uni10_double64* Mij_ori, int M, int N, uni10_double64* Q, uni10_double64* D, uni10_double64* R){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixLDQ(uni10_double64* Mij_ori, int M, int N, uni10_double64* L, uni10_double64* D, uni10_double64* Q){
      
      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixInv(uni10_double64* A, int N){

      uni10_error_msg(true, "%s", "Developping!!!\n");
      
    }

    uni10_double64 matrixDet(uni10_double64* A, int N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void setIdentity(uni10_double64* elem, size_t M, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void trimatrixEigH(uni10_double64* D, uni10_double64* E, int N, uni10_double64* z, int LDZ){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void eigSyDecompose(const uni10_double64* a, const uni10_int mp, const uni10_int nq, const uni10_int* desca, 
        uni10_int N, uni10_double64* w, uni10_double64* z, uni10_int* descz){
      
      uni10_int izero = 0, ione = 1;
      memcpy(z, a, N*N*sizeof(uni10_double64));
      
      uni10_int NPROW  = env_variables.get_info().nprow;
      uni10_int NPCOL  = env_variables.get_info().npcol;
      uni10_int NPROCS = env_variables.get_info().nprocs_mpi;
      uni10_int MYROW  = env_variables.get_info().myrow;

      //uni10_int nb = desca[5];
      //uni10_int nn = std::max(std::max(N, nb), 2, &NPROW);
      //uni10_int np = numroc(&nn, &nb, &izero, &izero, &NPROW);
      //uni10_int nq = numroc(&nn, &nb, &izero, &izero, &NPCOL);
      //uni10_int nrc = numroc(&N, &nb, &MYROW, &izero, &NPROCS);
      //uni10_int ldc = nrc;
      //uni10_int grmen = 2 * N - 2;
      //uni10_int lwmin = 5 * N + N * ldc + std::max(sizemqrleft, qrmen) + 1;


    }


  } /* namespace uni10_linalg */

} /* namespace uni10 */

