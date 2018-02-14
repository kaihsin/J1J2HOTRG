#ifndef __UNI10_SYS_INFO_CUSOLVER_GPU_H__
#define __UNI10_SYS_INFO_CUSOLVER_GPU_H__

#if defined(LINUX)
#include <sys/sysinfo.h>
#elif defined(OSX)
#include <sys/vmmeter.h>
#endif

#include <sys/types.h>
#include <sys/sysctl.h>

#include <unistd.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"


namespace uni10{

  enum uni10_runtime_type{
    only_cpu   =   0,
    hybrid     =   1,
    only_gpu   =   2
  };

  struct sysinfo_gpu{

    // Developping 
    sysinfo_gpu(): device(0), status(0){


      
    }

    void default_init(){

      status = 1;
      device = 0;

      cudaGetDeviceProperties(&devProp, device);

      globalMem          = devProp.totalGlobalMem;
      maxGridSize_x      = devProp.maxGridSize[0];
      maxGridSize_y      = devProp.maxGridSize[1];
      maxGridSize_z      = devProp.maxGridSize[2];
      maxThreadsPerBlock = devProp.maxThreadsPerBlock;
      maxThreadsDim_x    = devProp.maxThreadsDim[0];
      maxThreadsDim_y    = devProp.maxThreadsDim[1];
      maxThreadsDim_z    = devProp.maxThreadsDim[2];
      regsPerBlock       = devProp.regsPerBlock;

      // May be modified by the paramters in uni10rc.
      runtime_type       = only_gpu;
      threadsPerBlock_x  = 32;
      threadsPerBlock_y  = 32;

      checkCudaErrors(cudaMemGetInfo(&free_byte, &total_byte));

    }


    void clear(){


    };

    void print_sys_info() const{

      // Get device properties
      printf("\nCUDA Device #%d\n", device);
      cudaDeviceProp devProp;
      cudaGetDeviceProperties(&devProp, device);
      printf("Major revision number:         %d\n",  devProp.major);
      printf("Minor revision number:         %d\n",  devProp.minor);
      printf("Name:                          %s\n",  devProp.name);
      printf("Total global memory:           %ld\n",  devProp.totalGlobalMem);
      printf("Total shared memory per block: %ld\n",  devProp.sharedMemPerBlock);
      printf("Total registers per block:     %d\n",  devProp.regsPerBlock);
      printf("Warp size:                     %d\n",  devProp.warpSize);
      printf("Maximum memory pitch:          %ld\n",  devProp.memPitch);
      printf("Maximum threads per block:     %d\n",  devProp.maxThreadsPerBlock);
      for (int j = 0; j < 3; ++j)
        printf("Maximum dimension %d of block:  %d\n", j, devProp.maxThreadsDim[j]);
      for (int j = 0; j < 3; ++j)
        printf("Maximum dimension %d of grid:   %d\n", j, devProp.maxGridSize[j]);
      printf("Clock rate:                    %d\n",  devProp.clockRate);
      printf("Total constant memory:         %ld\n",  devProp.totalConstMem);
      printf("Texture alignment:             %ld\n",  devProp.textureAlignment);
      printf("Concurrent copy and execution: %s\n",  (devProp.deviceOverlap ? "Yes" : "No"));
      printf("Number of multiprocessors:     %d\n",  devProp.multiProcessorCount);
      printf("Kernel execution timeout:      %s\n",  (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
      printf("Compute:                       %d\n",  (devProp.computeMode));

    };

    uni10_int32 device;
    cudaDeviceProp devProp;

    uni10_uint64 globalMem, maxGridSize_x, maxGridSize_y, maxGridSize_z, maxThreadsPerBlock, maxThreadsDim_x, maxThreadsDim_y, maxThreadsDim_z, regsPerBlock;
    uni10_uint64 free_byte, total_byte;
    // The variable to handle the context of cublas.
    uni10_runtime_type runtime_type;

    uni10_int threadsPerBlock_x, threadsPerBlock_y;

    uni10_int status;  // status = 0. Has not initialized.

  };

}

#endif
