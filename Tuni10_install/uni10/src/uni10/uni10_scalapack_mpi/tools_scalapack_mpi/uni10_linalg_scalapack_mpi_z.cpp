#include <limits.h>

#include "uni10/uni10_error.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_tools_scalapack_mpi.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_linalg_scalapack_mpi_z.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_linalg_scalapack_mpi_dz.h"

#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_blacs_wrapper_mpi.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_scalapack_wrapper_mpi.h"

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(std::complex<double> a, std::complex<double>* X, int incx, std::complex<double>* Y, int incy, size_t N){   // Y = aX + Y

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    // Blas
    void vectorAdd(std::complex<double>* Y, std::complex<double>* X, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void vectorSub(std::complex<double>* Y, std::complex<double>* X, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void vectorMul(std::complex<double>* Y, std::complex<double>* X, size_t N){ 

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void vectorScal(std::complex<double> a, std::complex<double>* X, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void vectorExp(std::complex<double> a, std::complex<double>* X, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    std::complex<double> vectorSum(std::complex<double>* X, size_t N, int inc){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    double vectorNorm(std::complex<double>* X, size_t N, int inc){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixDot(std::complex<double>* A, std::complex<double>* B, int M, int N, int K, std::complex<double>* C){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void diagRowMul(std::complex<double>* mat, std::complex<double>* diag, size_t M, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void diagColMul(std::complex<double> *mat, std::complex<double>* diag, size_t M, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void setTranspose(std::complex<double>* A, size_t M, size_t N, std::complex<double>* AT){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void setTranspose(std::complex<double>* A, size_t M, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void setDagger(std::complex<double>* A, size_t M, size_t N, std::complex<double> *AT){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void setDagger(std::complex<double>* A, size_t M, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void setConjugate(std::complex<double> *A, size_t N, std::complex<double> *A_conj){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }


    void setConjugate(std::complex<double> *A, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    //LAPACK
    //
    void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, std::complex<double>* S_ori, std::complex<double>* vT){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixSDD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, std::complex<double>* S_ori, std::complex<double>* vT){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixQR(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* R){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixRQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* R, std::complex<double>* Q){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixLQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* L, std::complex<double>* Q){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixQL(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* L){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixQDR(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* D, std::complex<double>* R){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }


    void matrixLDQ(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* L, std::complex<double>* D, std::complex<double>* Q){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixQDRCPIVOT(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* Q, std::complex<double>* D, std::complex<double>* R){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }


    void matrixInv(std::complex<double>* A, int N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    std::complex<double> matrixDet(std::complex<double>* A, int N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void eigDecompose(std::complex<double>* Kij, int N, std::complex<double>* Eig, std::complex<double>* EigVec){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void setIdentity(std::complex<double>* elem, size_t M, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }


  } /* namespace uni10_linalg */

}
