#ifndef __UNI10_HIGH_RANK_LINALG_CONJ_H__
#define __UNI10_HIGH_RANK_LINALG_CONJ_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_conj.h"

namespace uni10{

  template<typename uni10_type>
    UniTensor<uni10_type> conj( const UniTensor<uni10_type>& T );

  template<typename uni10_type>
    UniTensor<uni10_type> conj( const UniTensor<uni10_type>& T ){

      UniTensor<uni10_type> Tout;
      conj(Tout, T, INPLACE);
      return Tout;

    }

};

#endif
