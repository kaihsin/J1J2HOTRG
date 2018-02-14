#include "uni10/uni10_sys_info/uni10_cusolver_gpu/uni10_cusolver_gpu_Dconst.cuh"
#include "uni10/uni10_env_info/uni10_cusolver_gpu/uni10_cusolver_gpu_initDconst.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/cuda_kernel_funcs/uni10_kernel_gpu.h"

namespace uni10{

  // uni10_double64
  __global__ void _vectorMul_kernel(uni10_double64* Y, uni10_double64* X, int N){

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < N)
      Y[i] *= X[i];

  }

  void vectorMul_kernel(uni10_double64* Y, uni10_double64* X, uni10_uint64 N){

    int64_t left        = N;
    uni10_uint64 offset = 0;
    uni10_int chunk;

    while(left > 0){

      if(left > INT_MAX)
        chunk = INT_MAX;
      else
        chunk = left;
      int64_t BLKPERGRID = (chunk + THREADSPERBLK_X - 1) / THREADSPERBLK_X;
      _vectorMul_kernel<<< BLKPERGRID, THREADSPERBLK_X>>>(Y + offset, X + offset, chunk);
      offset += chunk;
      left -= INT_MAX;

    }

  }

  __global__ void _vectorSum_kernel(uni10_double64 a, uni10_double64* X, uni10_uint64 N){

  }

  void vectorSum_kernel(uni10_double64 a, uni10_double64* X, uni10_uint64 N){

  }

  __global__ void _uni10_setDiag_kernel(uni10_double64* ori_elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N){


  }

  void uni10_setDiag_kernel(uni10_double64* ori_elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N){

  }

  __global__ void _uni10_getDiag_kernel(uni10_double64* ori_elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N){

  }

  void uni10_getDiag_kernel(uni10_double64* ori_elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N){

  }

  __global__ void _uni10_getUpTri(uni10_double64* ori_elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n){

  }

  void uni10_getUpTri(uni10_double64* ori_elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n){

  }

  __global__ void _uni10_getDnTri(uni10_double64* ori_elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n){

  }

  void uni10_getDnTri(uni10_double64* ori_elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n){

  }

  // uni10_complex128
  __global__ void _vectorMul_kernel(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N){

  }

  void vectorMul_kernel(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N){

  }

  __global__ void _vectorSum_kernel(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N){

  }

  void vectorSum_kernel(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N){

  }

  __global__ void _uni10_setDiag_kernel(uni10_complex128* ori_elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N){

  }

  void uni10_setDiag_kernel(uni10_complex128* ori_elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N){

  }

  __global__ void _uni10_getDiag_kernel(uni10_complex128* ori_elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N){

  }

  void uni10_getDiag_kernel(uni10_complex128* ori_elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N){

  }

  __global__ void _uni10_getUpTri(uni10_complex128* ori_elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n){

  }

  void uni10_getUpTri(uni10_complex128* ori_elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n){

  }

  __global__ void _uni10_getDnTri(uni10_complex128* ori_elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n){

  }

  void uni10_getDnTri(uni10_complex128* ori_elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n){

  }

  // Auxiliary double64-complex128 || complex128-double64
  __global__ void _uni10_elem_cast_kernel(uni10_complex128* new_elem, uni10_double64* raw_elem, uni10_uint64 elemNum){

  }

  void uni10_elem_cast_kernel(uni10_complex128* new_elem, uni10_double64* raw_elem, uni10_uint64 elemNum){

  }

  __global__ void _uni10_elem_cast_kernel(uni10_double64* new_elem, uni10_complex128* raw_elem, uni10_uint64 elemNum){

  }

  void uni10_elem_cast_kernel(uni10_double64* new_elem, uni10_complex128* raw_elem, uni10_uint64 elemNum){

  }

}
