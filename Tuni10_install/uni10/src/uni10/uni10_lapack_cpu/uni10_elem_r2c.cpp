#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  void Real2Complex(const uni10_elem_double64* d, uni10_elem_complex128* z){

    uni10_elem_cast(z->__elem, d->__elem, d->__elemNum);

  }

  void Real2Complex(const uni10_elem_complex128* _z, uni10_elem_complex128* z){

    uni10_error_msg(true, "%s", "Developing.");

  }
  
}
