#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  void setTranspose(const uni10_elem_double64* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* AT){

    uni10_linalg::setTranspose(A->__elem, *M, *N, AT->__elem);

  }

  void setTranspose(uni10_elem_double64* A, uni10_uint64* M, uni10_uint64* N){
    
    uni10_linalg::setTranspose(A->__elem, *M, *N);

  }

  void setTranspose(const uni10_elem_complex128* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* AT){

    uni10_linalg::setTranspose(A->__elem, *M, *N, AT->__elem);

  }

  void setTranspose(uni10_elem_complex128* A, uni10_uint64* M, uni10_uint64* N){

    uni10_linalg::setTranspose(A->__elem, *M, *N);

  }

}
