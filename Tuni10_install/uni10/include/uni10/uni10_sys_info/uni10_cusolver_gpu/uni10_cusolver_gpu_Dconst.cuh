#ifndef __UNI10_CUSOLVER_GPU_DCONST_CUH__
#define __UNI10_CUSOLVER_GPU_DCONST_CUH__

namespace uni10{

  __device__ __constant__ int MAXTHREADSPERBLOCK;
  __device__ __constant__ int MAXTHREADSDIM_X;
  __device__ __constant__ int MAXTHREADSDIM_Y;
  __device__ __constant__ int MAXTHREADSDIM_Z;
  __device__ __constant__ int MAXGRIDSIZE_X;
  __device__ __constant__ int MAXGRIDSIZE_Y;
  __device__ __constant__ int MAXGRIDSIZE_Z;

}; //End of uni10 namespace

#endif
