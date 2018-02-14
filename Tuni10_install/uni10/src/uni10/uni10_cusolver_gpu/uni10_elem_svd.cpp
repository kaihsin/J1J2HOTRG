#include "uni10/uni10_cusolver_gpu/uni10_elem_linalg_cusolver_gpu.h"

namespace uni10{

  void matrixSVD(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* U, uni10_elem_double64* S, uni10_elem_double64* vT){

    uni10_double64* U_elem  = (U  == NULL) ? NULL : U->__elem;
    uni10_double64* vT_elem = (vT == NULL) ? NULL : vT->__elem;

    if(!*isMdiag)
      uni10_linalg::matrixSVD(Mij_ori->__elem, *M, *N, U_elem, S->__elem, vT_elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

  void matrixSVD(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* U, uni10_elem_complex128* S, uni10_elem_complex128* vT){

    uni10_complex128* U_elem  = (U  == NULL) ? NULL : U->__elem;
    uni10_complex128* vT_elem = (vT == NULL) ? NULL : vT->__elem;

    if(!*isMdiag)
      uni10_linalg::matrixSVD(Mij_ori->__elem, *M, *N, U_elem, S->__elem, vT_elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

}
