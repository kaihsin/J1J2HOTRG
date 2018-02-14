#ifndef __UNI10_HIGH_RANK_LINALG_INPLACE_DAGGER_H__
#define __UNI10_HIGH_RANK_LINALG_INPLACE_DAGGER_H__

#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{

  template<typename uni10_type>
    void dagger( UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& Tin, UNI10_INPLACE on);

  template<typename uni10_type>
    void dagger( UniTensor<uni10_type>& T, UNI10_INPLACE on);

  template<typename uni10_type>
    void dagger( UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& Tin, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg(!(*Tin.status & UniTensor<uni10_type>::GET_HAVEBOND()), 
          "%s", "There is no bond in the tensor(scalar) to perform transposition.");

      uni10_uint64 bondNum = Tin.bonds->size();
      std::vector<uni10_int32> rsp_outin(bondNum);
      uni10_int32 rbondNum = 0;
      for(uni10_uint64 b = 0; b < bondNum; b++)
        if((*Tin.bonds)[b].type() == BD_IN)
          rbondNum++;
        else
          break;
      uni10_uint64 cbondNum = bondNum - rbondNum;
      for(uni10_uint64 b = 0; b < bondNum; b++)
        if(b < cbondNum)
          rsp_outin[b] = rbondNum + b;
        else
          rsp_outin[b] = b - cbondNum;
      std::vector<uni10_int32> outLabels(bondNum, 0);
      std::vector<Bond> outBonds;
      for(uni10_uint64 b = 0; b < Tin.bonds->size(); b++){
        outBonds.push_back((*Tin.bonds)[rsp_outin[b]]);
        outLabels[b] = (*Tin.labels)[rsp_outin[b]];
      }
      for(uni10_uint64 b = 0; b < bondNum; b++){
        if(b < cbondNum)
          outBonds[b].type_enforce() = BD_IN;
        else
          outBonds[b].type_enforce() = BD_OUT;
      }

      Tout.assign(outBonds);
      Tout.setName(*Tin.name);
      Tout.setLabel(outLabels);

      if((*Tin.status) & UniTensor<uni10_type>::GET_HAVEELEM())
        tensor_tools::dagger(Tout.paras, Tin.paras, Tin.style);

      *Tout.status |= UniTensor<uni10_type>::GET_HAVEELEM();

    }

  template<typename uni10_type>
    void dagger( UniTensor<uni10_type>& T, UNI10_INPLACE on){

      UniTensor<uni10_type> Tt;
      dagger(Tt, T, on);
      T = Tt;

    }

};

#endif
