#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"

namespace uni10{

  void vectorMul(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorMul(Y->__elem, X->__elem, *N);

  }

  void vectorMul(uni10_elem_complex128* Y, const uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorMul(Y->__elem, X->__elem, *N);

  }

  void vectorMul(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N){

    uni10_linalg::vectorMul(Y->__elem, X->__elem, *N);

  }  

}
