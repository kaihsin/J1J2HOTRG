#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  uni10_double64 matrixDet(const uni10_elem_double64* A, const uni10_uint64* N, uni10_const_bool* isMdiag){

    uni10_double64 det = 0.;
    if(!*isMdiag)
      det = uni10_linalg::matrixDet(A->__elem, *N);
    else
      uni10_error_msg(true, "%s", "Developping!!!");
    
    return det;  

  }

  uni10_complex128 matrixDet(const uni10_elem_complex128* A, const uni10_uint64* N, uni10_const_bool* isMdiag){

    uni10_complex128 det = 0.;
    if(!*isMdiag)
      det = uni10_linalg::matrixDet(A->__elem, *N);
    else
      uni10_error_msg(true, "%s", "Developping!!!");
    
    return det;  

  }

}
