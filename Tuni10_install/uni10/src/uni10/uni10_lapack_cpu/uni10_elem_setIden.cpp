#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  void setIdentity(uni10_elem_double64* A, const uni10_bool* is_diag, const uni10_uint64* M, const uni10_uint64* N){

    if(*is_diag == true){
      
      for(uni10_uint64 i = 0; i < A->__elemNum; i++)
        A->__elem[i] = 1.;

    }
    else{
      
      uni10_linalg::setIdentity(A->__elem, *M, *N);

    }

  }

  void setIdentity(uni10_elem_complex128* A, const uni10_bool* is_diag, const uni10_uint64* M, const uni10_uint64* N){

    if(*is_diag == true){
      
      for(uni10_uint64 i = 0; i < A->__elemNum; i++)
        A->__elem[i] = 1.;

    }
    else{
      
      uni10_linalg::setIdentity(A->__elem, *M, *N);

    }

  }

}
