#ifndef __UNI10_TOOLS_CUSOLVER_GPU_H__
#define __UNI10_TOOLS_CUSOLVER_GPU_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_env_info/uni10_cusolver_gpu/uni10_env_info_cusolver_gpu.h"

namespace uni10{

#ifdef __UNI10_ENV_CUSOLVER_GPU_H__
#define RUNTIMETYPE env_variables.get_info().runtime_type
#define UNI10_ONGPU env_variables.default_ongpu_flag()
#endif

  // Default values are going to be discarded.

  void* uni10_elem_alloc(uni10_uint64 memsize, uni10_bool& __ongpu);

  void* uni10_elem_copy(void* des, const void* src, uni10_uint64 memsize, uni10_bool des_ongpu = false, uni10_bool src_ongpu = false);

  void uni10_elem_free(void* ptr, uni10_uint64 memsize, uni10_bool __ongpu);

  void uni10_elemBzero(void* ptr, uni10_uint64 memsize, uni10_bool __ongpu);

  //For double ptr.
  //
  void uni10_setDiag(uni10_double64* ori_elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N, uni10_bool ori_ongpu = false, uni10_bool diag_ongpu = false);

  void uni10_getDiag(uni10_double64* ori_elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N, uni10_bool ori_ongpu = false, uni10_bool diag_ongpu = false);

  void uni10_print_elem_i(const uni10_double64& elem_i);

  void uni10_getUpTri(uni10_double64* ori_elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n, uni10_bool ori_ongpu = false, uni10_bool tri_ongpu = false); // elem -> tri_elem

  void uni10_getDnTri(uni10_double64* ori_elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n, uni10_bool ori_ongpu = false, uni10_bool tri_ongpu = false);
  //For complex ptr.
  //
  void uni10_setDiag(uni10_complex128* ori_elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N, uni10_bool ori_ongpu = false, uni10_bool diag_ongpu = false);

  void uni10_getDiag(uni10_complex128* ori_elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N, uni10_bool ori_ongpu = false, uni10_bool diag_ongpu = false);

  void uni10_print_elem_i(const uni10_complex128& elem_i);
  
  void uni10_getUpTri(uni10_complex128* ori_elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n, uni10_bool ori_ongpu = false, uni10_bool tri_ongpu = false); // elem -> tri_elem

  void uni10_getDnTri(uni10_complex128* ori_elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n, uni10_bool ori_ongpu = false, uni10_bool tri_ongpu = false);
  // Convert
  void uni10_elem_cast(uni10_complex128* des, const uni10_double64* src, uni10_uint64 N, uni10_bool des_ongpu = false, uni10_bool src_ongpu = false);

  void uni10_elem_cast(uni10_double64 *des, const uni10_complex128 *src, uni10_uint64 N, uni10_bool des_ongpu = false, uni10_bool src_ongpu = false);

  void ToReal(uni10_double64& M_i, uni10_double64 val);

  void ToReal(uni10_complex128& M_i, uni10_double64 val);

  void ToComplex(uni10_double64& M_i, uni10_double64 val);

  void ToComplex(uni10_complex128& M_i, uni10_double64 val);

  void shrinkWithoutFree(uni10_uint64 memsize);

}/* namespace uni10 */
#endif /* UNI10_TOOLS_H */
