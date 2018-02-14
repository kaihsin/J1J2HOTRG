#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"

namespace uni10{

  void setConjugate(const uni10_elem_double64* A, const uni10_uint64* N, uni10_elem_double64* A_conj){

    uni10_elem_copy(A_conj->__elem, A->__elem, *N *sizeof(uni10_elem_double64));

  }

  void setConjugate(uni10_elem_double64* A, uni10_uint64* N){

    //For function overload. Nothing to do.
    //
  }

  void setConjugate(const uni10_elem_complex128* A, const uni10_uint64* N, uni10_elem_complex128* A_conj){

    uni10_linalg::setConjugate(A->__elem, *N, A_conj->__elem);

  }

  void setConjugate(uni10_elem_complex128* A, uni10_uint64* N){

    uni10_linalg::setConjugate(A->__elem, *N);

  }

}
