#ifndef __UNI10_ELEM_RNG_H__
#define __UNI10_ELEM_RNG_H__

#if defined(UNI_LAPACK) && defined(UNI_CPU)
#include "uni10/uni10_lapack_cpu/uni10_elem_rng_lapack_cpu.h"
#endif

#if defined(UNI_CUSOLVER) && defined(UNI_GPU)
#include "uni10/uni10_cusolver_gpu/uni10_elem_rng_cusolver_gpu.h"
#endif

#if defined(UNI_SCALAPACK) && defined(UNI_MPI)
#include "uni10/uni10_scalapack_mpi/uni10_elem_rng_scalapack_mpi.h"
#endif

#endif
