#ifndef __UNI10_HIGH_RANK_LINALG_DAGGER_H__
#define __UNI10_HIGH_RANK_LINALG_DAGGER_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_dagger.h"

namespace uni10{

  template<typename uni10_type>
    UniTensor<uni10_type> dagger( const UniTensor<uni10_type>& T );

  template<typename uni10_type>
    UniTensor<uni10_type> dagger( const UniTensor<uni10_type>& T ){

      UniTensor<uni10_type> Tout;
      dagger(Tout, T, INPLACE);
      return Tout;

    }

};

#endif
