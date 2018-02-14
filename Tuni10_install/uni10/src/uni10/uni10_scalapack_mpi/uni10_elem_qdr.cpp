#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"

namespace uni10{

  void matrixQDR(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* D, uni10_elem_double64* R){
    
    if(!*isMdiag)
      uni10_linalg::matrixQDR(Mij_ori->__elem, *M, *N, Q->__elem, D->__elem, R->__elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

  void matrixQDR(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q,uni10_elem_complex128* D, uni10_elem_complex128* R){

    if(!*isMdiag)
      uni10_linalg::matrixQDR(Mij_ori->__elem, *M, *N, Q->__elem, D->__elem, R->__elem);
    else
      uni10_error_msg(true, "%s", "Developping!!!");

  }

}
