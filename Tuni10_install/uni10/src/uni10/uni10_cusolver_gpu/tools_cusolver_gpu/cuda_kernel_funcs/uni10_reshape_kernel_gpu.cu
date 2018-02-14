#include "uni10/uni10_sys_info/uni10_cusolver_gpu/uni10_cusolver_gpu_Dconst.cuh"
#include "uni10/uni10_env_info/uni10_cusolver_gpu/uni10_cusolver_gpu_initDconst.h"
#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/cuda_kernel_funcs/uni10_kernel_gpu.h"

namespace uni10{
  
  __global__ void _uni10_reshape_kernel(double* oldElem, int bondNum, size_t elemNum, size_t* offset, double* newElem){

    size_t oldIdx = blockIdx.y * MAXGRIDSIZE_X * MAXTHREADSPERBLOCK +  blockIdx.x * blockDim.x + threadIdx.x;
    size_t idx = oldIdx;
    size_t newIdx = 0;

    if(idx < elemNum){
      for(int i = 0; i < bondNum; i++){
        newIdx += (idx/offset[i]) * offset[bondNum + i];
        idx = idx % offset[i];
      }
      newElem[newIdx] = oldElem[oldIdx];
    }

  }

  void uni10_reshape_kernel(double* oldElem, int bondNum, size_t elemNum, size_t* offset, double* newElem){

    size_t* D_offset;
    checkCudaErrors(cudaMalloc((void**)&D_offset, 2 * sizeof(size_t) * bondNum));
    checkCudaErrors(cudaMemcpy(D_offset, offset, 2 * sizeof(size_t) * bondNum, cudaMemcpyHostToDevice));
    size_t blockNum = (elemNum + MAXTHREADSPERBLOCK_H - 1) / MAXTHREADSPERBLOCK_H;
    dim3 gridSize(blockNum % MAXGRIDSIZE_X_H, (blockNum + MAXGRIDSIZE_X_H - 1) / MAXGRIDSIZE_X_H);
    _uni10_reshape_kernel<<<gridSize, MAXTHREADSPERBLOCK_H>>>(oldElem, bondNum, elemNum, D_offset, newElem);

  }

  __global__ void _uni10_reshape_kernel(std::complex<double>* oldElem, int bondNum, size_t elemNum, size_t* offset, std::complex<double>* newElem){

    size_t oldIdx = blockIdx.y * MAXGRIDSIZE_X * MAXTHREADSPERBLOCK +  blockIdx.x * blockDim.x + threadIdx.x;
    size_t idx = oldIdx;
    size_t newIdx = 0;

    if(idx < elemNum){
      for(int i = 0; i < bondNum; i++){
        newIdx += (idx/offset[i]) * offset[bondNum + i];
        idx = idx % offset[i];
      }
      newElem[newIdx] = oldElem[oldIdx];
    }

  }

  void uni10_reshape_kernel(std::complex<double>* oldElem, int bondNum, size_t elemNum, size_t* offset, std::complex<double>* newElem){

    size_t* D_offset;
    checkCudaErrors(cudaMalloc((void**)&D_offset, 2 * sizeof(size_t) * bondNum));
    checkCudaErrors(cudaMemcpy(D_offset, offset, 2 * sizeof(size_t) * bondNum, cudaMemcpyHostToDevice));
    size_t blockNum = (elemNum + MAXTHREADSPERBLOCK_H - 1) / MAXTHREADSPERBLOCK_H;
    dim3 gridSize(blockNum % MAXGRIDSIZE_X_H, (blockNum + MAXGRIDSIZE_X_H - 1) / MAXGRIDSIZE_X_H);
    _uni10_reshape_kernel<<<gridSize, MAXTHREADSPERBLOCK_H>>>(oldElem, bondNum, elemNum, D_offset, newElem);

  }

}
