#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"

namespace uni10{

  void setDiag(uni10_elem_double64* _elem, const uni10_elem_double64* diag_elem, const uni10_uint64* M, const uni10_uint64* N ){

    uni10_setDiag(_elem->__elem, diag_elem->__elem, *M, *N, std::min(*M, *N));

  }

  void setDiag(uni10_elem_complex128* _elem, const uni10_elem_complex128* diag_elem, const uni10_uint64* M, const uni10_uint64* N ){

    uni10_setDiag(_elem->__elem, diag_elem->__elem, *M, *N, std::min(*M, *N));

  }

}
