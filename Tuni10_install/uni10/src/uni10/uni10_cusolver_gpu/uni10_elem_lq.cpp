#include "uni10/uni10_cusolver_gpu/uni10_elem_linalg_cusolver_gpu.h"

namespace uni10{

  void matrixLQ(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
       uni10_elem_double64* L, uni10_elem_double64* Q){

    if(!*isMdiag)
      uni10_linalg::matrixLQ(Mij_ori->__elem, *M, *N, L->__elem, Q->__elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

  void matrixLQ(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
       uni10_elem_complex128* L, uni10_elem_complex128* Q){

    if(!*isMdiag)
      uni10_linalg::matrixLQ(Mij_ori->__elem, *M, *N, L->__elem, Q->__elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

}
