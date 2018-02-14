#include <limits.h>

#include "uni10/uni10_error.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_tools_scalapack_mpi.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_linalg_scalapack_mpi_dz.h"

#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_blacs_wrapper_mpi.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_scalapack_wrapper_mpi.h"

namespace uni10{

  namespace uni10_linalg{

    void vectorAdd(double a, double* X, int incx, std::complex<double>* Y, int incy, size_t N){   // Y = aX + Y
      
      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void vectorAdd(std::complex<double>* Y, double* X, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void vectorSub(std::complex<double>* Y, double* X, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void vectorMul(std::complex<double>* Y, double* X, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void vectorScal(double a, std::complex<double>* X, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void vectorExp(double a, std::complex<double>* X, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixDot(double* A, std::complex<double>* B, int M, int N, int K, std::complex<double>* C){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixDot(std::complex<double>* A, double* B, int M, int N, int K, std::complex<double>* C){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void diagRowMul(std::complex<double>* mat, double* diag, size_t M, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void diagColMul(std::complex<double>* mat, double* diag, size_t M, size_t N){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void eigDecompose(double* Kij_ori, int N, std::complex<double>* Eig, std::complex<double>* EigVec){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void eigSyDecompose(std::complex<double>* Kij, int N, double* Eig, std::complex<double>* EigVec){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixSVD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, double *S, std::complex<double>* vT){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

    void matrixSDD(std::complex<double>* Mij_ori, int M, int N, std::complex<double>* U, double *S, std::complex<double>* vT){

      uni10_error_msg(true, "%s", "Developping!!!\n");

    }

  };/* namespace uni10_linalg */

};/* namespace uni10 */

