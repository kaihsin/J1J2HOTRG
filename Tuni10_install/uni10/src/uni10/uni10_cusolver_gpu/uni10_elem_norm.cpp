#include "uni10/uni10_cusolver_gpu/uni10_elem_linalg_cusolver_gpu.h"

namespace uni10{

  uni10_double64 vectorNorm(const uni10_elem_double64* X, const uni10_uint64* N, uni10_int32* inc){

    return uni10_linalg::vectorNorm(X->__elem, *N, *inc);

  }

  uni10_double64 vectorNorm(const uni10_elem_complex128* X, const uni10_uint64* N, uni10_int32* inc){

    return uni10_linalg::vectorNorm(X->__elem, *N, *inc);

  }
  
}
