#ifndef __UNI10_HIGH_RANK_LINALG_TRACE_H__
#define __UNI10_HIGH_RANK_LINALG_TRACE_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_permute.h"
#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{

  template<typename uni10_type>
    uni10_type trace(const UniTensor<uni10_type>& T);

  template<typename uni10_type>
    uni10_type trace(const UniTensor<uni10_type>& T){

      uni10_error_msg(!((*T.status) & T.HAVEELEM), "%s", "Cannot trace a tensor before setting its elements.");

      uni10_type trVal = 0.;
      if((*T.status) & T.HAVEBOND){
        typename std::map<Qnum, Block<uni10_type> >::const_iterator it = T.blocks->begin();
        for( ; it != T.blocks->end(); it++ ){
          uni10_error_msg(!(it->second.row_enforce() == it->second.col_enforce()), "%s", "Cannot trace a non-square block.");
          trVal += trace(it->second);
        }
        return trVal;
      }else
        return trVal;
    }  

};

#endif
