#include <chrono>

#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"
#include "uni10/uni10_scalapack_mpi/uni10_elem_rng_scalapack_mpi.h"

namespace uni10{

  void setUniformRand(uni10_elem_double64* A, const uni10_bool* is_diag, const uni10_uint64* M, const uni10_uint64* N,
      const uni10_double64* up, const uni10_double64* dn, const uni10_int64* _seed){

    uni10_uint64 seed = (*_seed < 0) ? std::chrono::system_clock::now().time_since_epoch().count() : (uni10_uint64)*_seed; 
    uni10_mt19937 generator(seed);
    uni10_uniform_real distribution(*up, *dn);

    uni10_uint64 elemNum = *is_diag ? std::min(*M, *N) : (*M)*(*N);

    for(uni10_uint64 i = 0; i < elemNum; i++)
      A->__elem[i] = distribution(generator);

  }

  void setUniformRand(uni10_elem_complex128* A, const uni10_bool* is_diag, const uni10_uint64* M, const uni10_uint64* N,
      const uni10_double64* up, const uni10_double64* dn, const uni10_int64* _seed){

    uni10_uint64 seed = (*_seed < 0) ? std::chrono::system_clock::now().time_since_epoch().count() : (uni10_uint64)*_seed; 
    uni10_mt19937 generator(seed);
    uni10_uniform_real distribution(*up, *dn);

    uni10_uint64 elemNum = *is_diag ? std::min(*M, *N) : (*M)*(*N);

    for(uni10_uint64 i = 0; i < elemNum; i++){
      A->__elem[i].real(distribution(generator));
      A->__elem[i].imag(distribution(generator));
    }

  }

}
