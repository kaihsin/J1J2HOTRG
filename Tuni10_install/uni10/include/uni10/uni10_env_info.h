#ifndef __UNI10_ENV_H__
#define __UNI10_ENV_H__ 

#if defined(UNI_CPU) && defined(UNI_LAPACK)
#include "uni10/uni10_env_info/uni10_lapack_cpu/uni10_env_info_lapack_cpu.h"
#endif

#if defined(UNI_GPU) && defined(UNI_CUSOLVER)
#include "uni10/uni10_env_info/uni10_cusolver_gpu/uni10_env_info_cusolver_gpu.h"
#endif

#if defined(UNI_MPI) && defined(UNI_SCALAPACK)
#include "uni10/uni10_env_info/uni10_scalapack_mpi/uni10_env_info_scalapack_mpi.h"
#endif

#endif
