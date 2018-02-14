#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"

namespace uni10{

  void matrixInv(const uni10_elem_double64* A, const uni10_uint64* N, uni10_const_bool* isMdiag){

    if(!*isMdiag)
      uni10_linalg::matrixInv(A->__elem, *N);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

  void matrixInv(const uni10_elem_complex128* A, const uni10_uint64* N, uni10_const_bool* isMdiag){

    if(!*isMdiag)
      uni10_linalg::matrixInv(A->__elem, *N);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

}
