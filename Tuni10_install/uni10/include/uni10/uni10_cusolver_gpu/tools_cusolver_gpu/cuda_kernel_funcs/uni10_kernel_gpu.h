#ifndef __UNI10_KERNEL_GPU_H__
#define __UNI10_KERNEL_GPU_H__

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_env_info.h"

namespace uni10{

#ifdef __UNI10_ENV_CUSOLVER_GPU_H__
#ifndef THREADSPERBLK_X
#define THREADSPERBLK_X env_variables.get_info().threadsPerBlock_x
#endif
#ifndef THREADSPERBLK_Y
#define THREADSPERBLK_Y env_variables.get_info().threadsPerBlock_y
#endif
#ifndef MAXGRIDSIZE_X_H
#define MAXGRIDSIZE_X_H env_variables.get_info().maxGridSize_x
#endif
#ifndef MAXGRIDSIZE_Y_H
#define MAXGRIDSIZE_Y_H env_variables.get_info().maxGridSize_y
#endif
#ifndef MAXGRIDSIZE_Z_H
#define MAXGRIDSIZE_Z_H env_variables.get_info().maxGridSize_z
#endif
#ifndef MAXTHREADSPERBLOCK_H
#define MAXTHREADSPERBLOCK_H env_variables.get_info().maxThreadsPerBlock
#endif
#endif

  // uni10_double64
  void vectorExp_kernel(uni10_double64 a, uni10_double64* X, uni10_uint64 N);

  //element-wise multiplication.
  void vectorMul_kernel(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);

  void vectorSum_kernel(uni10_double64 a, uni10_double64* X, uni10_uint64 N);

  void uni10_setDiag_kernel(uni10_double64* ori_elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_getDiag_kernel(uni10_double64* ori_elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_getUpTri(uni10_double64* ori_elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n); // elem -> tri_elem

  void uni10_getDnTri(uni10_double64* ori_elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n);

  // uni10_complex128
  void vectorExp_kernel(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);

  //element-wise multiplication.
  void vectorMul_kernel(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);

  void vectorSum_kernel(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);

  void uni10_setDiag_kernel(uni10_complex128* ori_elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_getDiag_kernel(uni10_complex128* ori_elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_getUpTri(uni10_complex128* ori_elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n); // elem -> tri_elem

  void uni10_getDnTri(uni10_complex128* ori_elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n);

  // Auxiliary double64-complex128 || complex128-double64
  void uni10_elem_cast_kernel(uni10_complex128* new_elem, uni10_double64* raw_elem, uni10_uint64 elemNum);

  void uni10_elem_cast_kernel(uni10_double64* new_elem, uni10_complex128* raw_elem, uni10_uint64 elemNum);

  void uni10_reshape_kernel(double* oldElem, int bondNum, size_t elemNum, size_t* offset, double* newElem);

  void uni10_reshape_kernel(std::complex<double>* oldElem, int bondNum, size_t elemNum, size_t* offset, std::complex<double>* newElem);

}

#endif
