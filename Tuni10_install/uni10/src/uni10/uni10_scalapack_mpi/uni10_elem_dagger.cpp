#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"

namespace uni10{

  void setDagger(const uni10_elem_double64* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* AT){

    uni10_linalg::setDagger(A->__elem, A->desc, *M, *N, AT->__elem, AT->desc);

  }

  void setDagger(uni10_elem_double64* A, uni10_uint64* M, uni10_uint64* N){

    uni10_elem_double64 AT(*N, *M, A->__elemNum != A->Rnum_ * A->Cnum_);
    uni10_linalg::setDagger(A->__elem, A->desc, *M, *N, AT.__elem, AT.desc);
    *A = AT;

  }

  void setDagger(const uni10_elem_complex128* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* AT){

    uni10_linalg::setDagger(A->__elem, *M, *N, AT->__elem);

  }

  void setDagger(uni10_elem_complex128* A, uni10_uint64* M, uni10_uint64* N){

    uni10_linalg::setDagger(A->__elem, *M, *N);

  }

}
