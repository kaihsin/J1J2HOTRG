#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  void vectorAdd(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorAdd(Y->__elem, X->__elem, *N);

  }

  void vectorAdd(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N){

    uni10_linalg::vectorAdd(Y->__elem, X->__elem, *N);

  }

  void vectorAdd(uni10_elem_complex128* Y, const uni10_elem_double64* X, const uni10_uint64* N){

    uni10_linalg::vectorAdd(Y->__elem, X->__elem, *N);

  }


}
