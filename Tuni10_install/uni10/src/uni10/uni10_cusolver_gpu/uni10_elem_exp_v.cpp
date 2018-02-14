#include "uni10/uni10_cusolver_gpu/uni10_elem_linalg_cusolver_gpu.h"

namespace uni10{

  void vectorExp(uni10_double64* a, uni10_elem_double64* X, uni10_uint64* N){

    uni10_linalg::vectorExp(*a, X->__elem, *N);

  }

  void vectorExp(uni10_complex128* a, uni10_elem_complex128* X, uni10_uint64* N){

    uni10_linalg::vectorExp(*a, X->__elem, *N);

  }
  
  void vectorExp(uni10_double64* a, uni10_elem_complex128* X, uni10_uint64* N){

    uni10_linalg::vectorExp(*a, X->__elem, *N);

  }

}
