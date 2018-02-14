#ifndef __UNI10_HIGH_RANK_LINALG_CONTRACT_H__
#define __UNI10_HIGH_RANK_LINALG_CONTRACT_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_contract.h"

namespace uni10{

  template<typename T>
    UniTensor<T> contract(UniTensor<uni10_double64>& Ta, UniTensor<T>& Tb, bool fast);

  template<typename T>
    UniTensor<uni10_complex128> contract(UniTensor<uni10_complex128>& Ta, UniTensor<T>& Tb, bool fast);

  template<typename T>
    UniTensor<T> contract(UniTensor<uni10_double64>& Ta, UniTensor<T>& Tb, bool fast){

      UniTensor<T> Tout;
      contract(Tout, Ta, Tb, fast, INPLACE);
      return Tout;

    }

  template<typename T>
    UniTensor<uni10_complex128> contract(UniTensor<uni10_complex128>& Ta, UniTensor<T>& Tb, bool fast){

      UniTensor<uni10_complex128> Tout;
      contract(Tout, Ta, Tb, fast, INPLACE);
      return Tout;

    }

};

#endif
