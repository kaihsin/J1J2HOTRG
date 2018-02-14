#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"

namespace uni10{

  uni10_double64 matrixTrace(const uni10_elem_double64* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N){
    
    uni10_uint64 min = std::min(*M, *N);
    uni10_int32 inc = 1;
    uni10_double64 res = 0.;
    if(!*isMdiag){
      for(uni10_uint64 i = 0; i < min; i++){
        res += Mij_ori->__elem[i*(*N)+i];
      }
    }
    else
      res = uni10_linalg::vectorSum(Mij_ori->__elem, min, inc);
    return res;

  }

  uni10_complex128 matrixTrace(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isMdiag, const uni10_uint64* M, const uni10_uint64* N){

    uni10_uint64 min = std::min(*M, *N);
    uni10_int32 inc = 1;
    uni10_complex128 res = 0.;
    if(!*isMdiag){
      for(uni10_uint64 i = 0; i < min; i++){
        res += Mij_ori->__elem[i*(*N)+i];
      }
    }
    else
      res = uni10_linalg::vectorSum(Mij_ori->__elem, min, inc);
    return res;
  }

}
