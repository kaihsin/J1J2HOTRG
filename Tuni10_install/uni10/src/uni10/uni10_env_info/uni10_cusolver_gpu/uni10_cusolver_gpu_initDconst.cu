#include "uni10/uni10_sys_info/uni10_cusolver_gpu/uni10_cusolver_gpu_Dconst.cuh"
#include "uni10/uni10_env_info/uni10_cusolver_gpu/uni10_cusolver_gpu_initDconst.h"

namespace uni10{

  void uni10_init_device_const(sysinfo_gpu& uni10_sys_info){

    uni10_int memsize = sizeof(uni10_int);
    cudaMemcpyToSymbol( MAXGRIDSIZE_X, &uni10_sys_info.maxGridSize_x, memsize);
    cudaMemcpyToSymbol( MAXGRIDSIZE_Y, &uni10_sys_info.maxGridSize_y, memsize);
    cudaMemcpyToSymbol( MAXGRIDSIZE_Z, &uni10_sys_info.maxGridSize_z, memsize);
    cudaMemcpyToSymbol( MAXTHREADSPERBLOCK, &uni10_sys_info.maxThreadsPerBlock, memsize);
    cudaMemcpyToSymbol( MAXTHREADSDIM_X, &uni10_sys_info.maxThreadsDim_x, memsize);
    cudaMemcpyToSymbol( MAXTHREADSDIM_Y, &uni10_sys_info.maxThreadsDim_y, memsize);
    cudaMemcpyToSymbol( MAXTHREADSDIM_Z, &uni10_sys_info.maxThreadsDim_z, memsize);

  }

  std::map<std::string, uni10_int> uni10_get_device_const(){
    
    uni10_int memsize = sizeof(uni10_int);

    std::map<std::string, uni10_int> dev_info;
    uni10_int para;
    cudaMemcpyFromSymbol( &para, MAXGRIDSIZE_X, memsize);
    dev_info["MAXGRIDSIZE_X"] = para;
    cudaMemcpyFromSymbol( &para, MAXGRIDSIZE_Y, memsize);
    dev_info["MAXGRIDSIZE_Y"] = para;
    cudaMemcpyFromSymbol( &para, MAXGRIDSIZE_Z, memsize);
    dev_info["MAXGRIDSIZE_Z"] = para;
    cudaMemcpyFromSymbol( &para, MAXTHREADSPERBLOCK, memsize);
    dev_info["MAXTHREADSPERBLOCK"] = para;
    cudaMemcpyFromSymbol( &para, MAXTHREADSDIM_X, memsize);
    dev_info["MAXTHREADSDIM_X"] = para;
    cudaMemcpyFromSymbol( &para, MAXTHREADSDIM_Y, memsize);
    dev_info["MAXTHREADSDIM_Y"] = para;
    cudaMemcpyFromSymbol( &para, MAXTHREADSDIM_Z, memsize);
    dev_info["MAXTHREADSDIM_Z"] = para;

    return dev_info;

  }

}; // End of uni10 namespace

