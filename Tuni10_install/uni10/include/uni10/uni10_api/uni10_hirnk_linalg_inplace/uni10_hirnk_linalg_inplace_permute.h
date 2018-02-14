#ifndef __UNI10_HIGH_RANK_LINALG_INPLACE_PERMUTE_H__
#define __UNI10_HIGH_RANK_LINALG_INPLACE_PERMUTE_H__

#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& T, const std::vector<uni10_int32>& newLabels, uni10_int32 rowBondNum, UNI10_INPLACE on);

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& T, uni10_int32* newLabels, uni10_int32 rowBondNum, UNI10_INPLACE on);

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& T, uni10_int32 rowBondNum, UNI10_INPLACE on);

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& T, const std::vector<uni10_int32>& newLabels, uni10_int32 rowBondNum, UNI10_INPLACE on);

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& T, uni10_int32* newLabels, uni10_int32 rowBondNum, UNI10_INPLACE on);

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& T, uni10_int32 rowBondNum, UNI10_INPLACE on);


  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& T, const std::vector<uni10_int32>& newLabels, uni10_int32 rowBondNum, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg((*T.status) & UniTensor<uni10_type>::GET_HAVEBOND == 0, "%s", "There is no bond in the tensor(scalar) to permute.");
      uni10_error_msg((T.labels->size() == newLabels.size()) == 0, "%s", "The size of the input new labels does not match for the number of bonds.");

      uni10_int32 bondNum = T.bonds->size();
      std::vector<uni10_int32> rsp_outin(bondNum);
      uni10_int32 cnt = 0;
      for(uni10_int32 i = 0; i < bondNum; i++)
        for(uni10_int32 j = 0; j < bondNum; j++)
          if((*T.labels)[i] == newLabels[j]){
            rsp_outin[j] = i;
            cnt++;
          }
      uni10_error_msg((cnt == newLabels.size()) == 0, "%s", "The input new labels do not 1-1 correspond to the labels of the tensor.");

      uni10_bool inorder = true;

      for(uni10_int32 i = 1; i < bondNum; i++)
        if(rsp_outin[i] != i){
          inorder = false;
          break;
        }
      if(inorder && (*T.RBondNum) == rowBondNum) {
        Tout = T;
        return ;
      } 
      else{
        std::vector<Bond> outBonds;
        for(uni10_int32 b = 0; b < T.bonds->size(); b++){
          outBonds.push_back((*T.bonds)[rsp_outin[b]]);
        }
        for(uni10_uint64 b = 0; b < T.bonds->size(); b++){
          if(b < (uni10_uint64)rowBondNum)
            outBonds[b].change(BD_IN);
          else
            outBonds[b].change(BD_OUT);
        }
        Tout.assign(outBonds);
        Tout.setName(*T.name);
        // ON CPU 
        // ON GPU Developping
        if((*T.status) & UniTensor<uni10_type>::GET_HAVEELEM())
          tensor_tools::permute(T.paras, T.style, rsp_outin, Tout.paras, inorder);

        *Tout.status |= UniTensor<uni10_type>::GET_HAVEELEM();
        Tout.setLabel(newLabels);

      }
    }

  template<typename uni10_type>
    void permute(  UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& T, uni10_int32* newLabels, uni10_int32 rowBondNum, UNI10_INPLACE on){

      std::vector<uni10_int32> _labels(newLabels, newLabels + T.bond().size());
      permute( Tout, T, _labels, rowBondNum, on);

    }

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& T, uni10_int32 rowBondNum, UNI10_INPLACE on){

      std::vector<uni10_int32> ori_labels = T.label();
      permute( Tout, T, ori_labels, rowBondNum, on);

    }

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& T, const std::vector<uni10_int32>& newLabels, uni10_int32 rowBondNum, UNI10_INPLACE on){
      UniTensor<uni10_type> Tout;
      permute( Tout, T, newLabels, rowBondNum, on);
      T = Tout;
    }

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& T, uni10_int32* newLabels, uni10_int32 rowBondNum, UNI10_INPLACE on){

      UniTensor<uni10_type> Tout;
      permute( Tout, T, newLabels, rowBondNum, on);
      T = Tout;

    }

  template<typename uni10_type>
    void permute( UniTensor<uni10_type>& T, uni10_int32 rowBondNum, UNI10_INPLACE on){

      UniTensor<uni10_type> Tout;
      permute( Tout, T, rowBondNum, on);
      T = Tout;

    }

};

#endif
