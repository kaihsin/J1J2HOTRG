#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"

namespace uni10{

  void matrixQL(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* L){

    if(!*isMdiag)
      uni10_linalg::matrixQL(Mij_ori->__elem, Mij_ori->blockrow, Mij_ori->blockcol, Mij_ori->desc, *M, *N, Q->__elem, Q->desc, L->__elem, L->desc);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

  void matrixQL(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* L){

    if(!*isMdiag)
      uni10_linalg::matrixQL(Mij_ori->__elem, *M, *N, Q->__elem, L->__elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

}
