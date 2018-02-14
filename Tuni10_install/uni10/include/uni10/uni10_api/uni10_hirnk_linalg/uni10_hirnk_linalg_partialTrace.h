#ifndef __UNI10_HIGH_RANK_LINALG_PARTIALTRACE_H__
#define __UNI10_HIGH_RANK_LINALG_PARTIALTRACE_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_permute.h"
#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{

  template<typename uni10_type>
    UniTensor<uni10_type> partialTrace(const UniTensor<uni10_type>& T, uni10_int32 la, uni10_int32 lb);

  template<typename uni10_type>
    UniTensor<uni10_type> partialTrace(const UniTensor<uni10_type>& T, uni10_int32 la, uni10_int32 lb){

      uni10_error_msg(!((*T.status) & T.HAVEELEM), "%s" ,"Cannot perform contraction of two tensors before setting their elements.");
      uni10_error_msg(!(T.bonds->size() > 2), "%s", "The number of bonds must larger than 2 for performing partialTrace.");

      uni10_int32 bondNum = T.bonds->size();
      std::vector<Bond> newBonds;
      std::vector<uni10_int32>newLabels(bondNum - 2, 0);
      std::vector<uni10_int32>rsp_labels(bondNum);
      uni10_int32 ia, ib;
      uni10_int32 enc = 0;

      for(uni10_uint64 l = 0; l < T.labels->size(); l++){
        if((*T.labels)[l] == la)
          ia = l;
        else if((*T.labels)[l] == lb)
          ib = l;
        else{
          newBonds.push_back((*T.bonds)[l]);
          newLabels[enc] = (*T.labels)[l];
          rsp_labels[enc] = (*T.labels)[l];
          enc++;
        }
      }

      uni10_error_msg(!(enc == newLabels.size()), "%s", "Cannot find the two bonds with the given two labels.");

      UniTensor<uni10_type> Tt(newBonds, newLabels);
      rsp_labels[bondNum - 2] = (*T.labels)[ia];
      rsp_labels[bondNum - 1] = (*T.labels)[ib];
      ia = bondNum - 2;
      ib = bondNum - 1;

      UniTensor<uni10_type> Tin;
      permute(Tin, T, rsp_labels, *Tt.RBondNum, INPLACE);

      tensor_tools::traceByRow(Tt.paras, Tin.paras, ia, ib, Tin.style);

      *Tt.status |= UniTensor<uni10_type>::GET_HAVEELEM();

      return Tt;

    }

};

#endif
