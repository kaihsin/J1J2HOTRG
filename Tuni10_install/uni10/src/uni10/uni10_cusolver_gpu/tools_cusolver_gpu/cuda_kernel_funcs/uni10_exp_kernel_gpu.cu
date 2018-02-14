#include "uni10/uni10_sys_info/uni10_cusolver_gpu/uni10_cusolver_gpu_Dconst.cuh"
#include "uni10/uni10_env_info/uni10_cusolver_gpu/uni10_cusolver_gpu_initDconst.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/cuda_kernel_funcs/uni10_kernel_gpu.h"

namespace uni10{
  
  __global__ void _vectorExp_kernel(uni10_double64 a, uni10_double64* X, uni10_uint64 N){

    uni10_uint64 idx = blockIdx.y * MAXGRIDSIZE_X * MAXTHREADSPERBLOCK + blockIdx.x * blockDim.x + threadIdx.x;

    if(idx < N)
      X[idx] = exp( a * X[idx]);

  }

  void vectorExp_kernel(uni10_double64 a, uni10_double64* X, uni10_uint64 N){

    uni10_uint64 blockNum = (N + MAXTHREADSPERBLOCK_H - 1) / MAXTHREADSPERBLOCK_H;
    dim3 gridSize(blockNum % MAXGRIDSIZE_X_H, (blockNum + MAXGRIDSIZE_X_H - 1) / MAXGRIDSIZE_X_H);
    _vectorExp_kernel<<<gridSize, MAXTHREADSPERBLOCK_H>>>(a, X, N);

  }

  __global__ void _vectorExp_kernel(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N){

    //uni10_uint64 idx = blockIdx.y * MAXGRIDSIZE_X * MAXTHREADSPERBLOCK + blockIdx.x * blockDim.x + threadIdx.x;

    //if(idx < N)
    //  X[idx] = exp( a * X[idx]);

  }

  void vectorExp_kernel(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N){

    uni10_uint64 blockNum = (N + MAXTHREADSPERBLOCK_H - 1) / MAXTHREADSPERBLOCK_H;
    dim3 gridSize(blockNum % MAXGRIDSIZE_X_H, (blockNum + MAXGRIDSIZE_X_H - 1) / MAXGRIDSIZE_X_H);
    _vectorExp_kernel<<<gridSize, MAXGRIDSIZE_X_H>>>(a, X, N);
    uni10_error_msg(true, "%s", "Developing");

  }

}
