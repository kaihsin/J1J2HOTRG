#ifndef __UNI10_HIGH_RANK_LINALG_TRANS_H__
#define __UNI10_HIGH_RANK_LINALG_TRANS_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_trans.h"

namespace uni10{

  template<typename uni10_type>
    UniTensor<uni10_type> transpose( const UniTensor<uni10_type>& T );

  template<typename uni10_type>
    UniTensor<uni10_type> transpose( const UniTensor<uni10_type>& T ){

      UniTensor<uni10_type> Tout;
      transpose(Tout, T, INPLACE);
      return Tout;

    }


};

#endif
