#include "uni10/uni10_cusolver_gpu/uni10_elem_linalg_cusolver_gpu.h"

namespace uni10{

  void matrixEigh(const uni10_elem_double64* Mij_ori, uni10_const_bool *isMdiag, const uni10_uint64* N, uni10_elem_double64* Eig, uni10_elem_double64* EigVec){

    if(!*isMdiag)
      uni10_linalg::eigSyDecompose(Mij_ori->__elem, *N, Eig->__elem, EigVec->__elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

  void matrixEigh(const uni10_elem_complex128* Mij_ori, uni10_const_bool *isMdiag, const uni10_uint64* N, uni10_elem_double64* Eig, uni10_elem_complex128* EigVec){

    if(!*isMdiag)
      uni10_linalg::eigSyDecompose(Mij_ori->__elem, *N, Eig->__elem, EigVec->__elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

  void matrixEigh(const uni10_elem_complex128* Mij_ori, uni10_const_bool *isMdiag, const uni10_uint64* N, uni10_elem_complex128* Eig, uni10_elem_complex128* EigVec){

    uni10_double64* tmp = (uni10_double64*)malloc((*N)*sizeof(uni10_double64));
    if(!*isMdiag){
      uni10_linalg::eigSyDecompose(Mij_ori->__elem, *N, tmp, EigVec->__elem);
      uni10_elem_cast(Eig->__elem, tmp, (*N));
    }
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

}
