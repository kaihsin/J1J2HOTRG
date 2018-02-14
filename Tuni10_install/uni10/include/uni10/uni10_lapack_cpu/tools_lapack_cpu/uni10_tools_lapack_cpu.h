#ifndef __UNI10_TOOLS_CPU_H__
#define __UNI10_TOOLS_CPU_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"

#ifdef UNI_MKL
#include "mkl.h"
#else
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_lapack_wrapper_cpu.h"
#endif

namespace uni10{

  void* uni10_elem_alloc(uni10_uint64 memsize);

  void* uni10_elem_copy(void* des, const void* src, uni10_uint64 memsize);

  void uni10_elem_free(void* ptr, uni10_uint64 memsize);

  void uni10_elemBzero(void* ptr, uni10_uint64 memsize);

  //For double ptr.
  //
  void uni10_setDiag(uni10_double64* elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_getDiag(uni10_double64* elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_print_elem_i(const uni10_double64& elem_i);

  void uni10_getUpTri(uni10_double64* elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n); // elem -> tri_elem

  void uni10_getDnTri(uni10_double64* elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n);
  //For complex ptr.
  //
  void uni10_setDiag(uni10_complex128* elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_getDiag(uni10_complex128* elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_print_elem_i(const uni10_complex128& elem_i);
  
  void uni10_getUpTri(uni10_complex128* elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n); // elem -> tri_elem

  void uni10_getDnTri(uni10_complex128* elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n);
  // Convert
  void uni10_elem_cast(uni10_complex128* des, const uni10_double64* src, uni10_uint64 N);

  void uni10_elem_cast(uni10_double64 *des, const uni10_complex128 *src, uni10_uint64 N);

  void ToReal(uni10_double64& M_i, uni10_double64 val);

  void ToReal(uni10_complex128& M_i, uni10_double64 val);

  void ToComplex(uni10_double64& M_i, uni10_double64 val);

  void ToComplex(uni10_complex128& M_i, uni10_double64 val);

  void shrinkWithoutFree(uni10_uint64 memsize);

}/* namespace uni10 */
#endif /* UNI10_TOOLS_H */
