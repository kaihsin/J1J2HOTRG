#include "uni10/uni10_cusolver_gpu/uni10_elem_linalg_cusolver_gpu.h"

namespace uni10{

  void trimatrixEigh(uni10_elem_double64* D, uni10_elem_double64* E, uni10_uint64* N, 
      uni10_elem_double64* z, uni10_uint64* LDZ){

    if(z==NULL && LDZ==NULL)
      uni10_linalg::trimatrixEigH(D->__elem, E->__elem, *N);
    else
      uni10_linalg::trimatrixEigH(D->__elem, E->__elem, *N, z->__elem, *LDZ);

  }

  void trimatrixEigh(uni10_elem_complex128* D, uni10_elem_complex128* E, uni10_uint64* N, 
      uni10_elem_complex128* z, uni10_uint64* LDZ){

    uni10_error_msg(true, "%s", "developping!!");

  }

}
