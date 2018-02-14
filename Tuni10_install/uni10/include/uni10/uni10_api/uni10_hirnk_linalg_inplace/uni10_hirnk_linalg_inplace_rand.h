#ifndef __UNI10_HIGH_RANK_LINALG_INPLACE_RAND_H__
#define __UNI10_HIGH_RANK_LINALG_INPLACE_RAND_H__

#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{

  template<typename uni10_type> 
    void randomize( UniTensor<uni10_type>& A );

  template<typename uni10_type> 
    void randomize( UniTensor<uni10_type>& A ){

      tensor_tools::randomize(A.paras, A.style);

    }

};

#endif
