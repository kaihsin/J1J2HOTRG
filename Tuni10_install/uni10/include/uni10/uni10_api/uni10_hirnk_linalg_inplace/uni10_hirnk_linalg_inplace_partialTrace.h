#ifndef __UNI10_HIGH_RANK_LINALG_INPLACE_PARTIALTRA_H__
#define __UNI10_HIGH_RANK_LINALG_INPLACE_PARTIALTRA_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_permute.h"
#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{

  template<typename uni10_type>
    void partialTrace(UniTensor<uni10_type>& T, uni10_int32 la, uni10_int32 lb, UNI10_INPLACE on);

  template<typename uni10_type>
    void partialTrace(UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& Tin, uni10_int32 la, uni10_int32 lb, UNI10_INPLACE on);

  template<typename uni10_type>
    void partialTrace(UniTensor<uni10_type>& T, uni10_int32 la, uni10_int32 lb, UNI10_INPLACE on){

      UniTensor<uni10_type> Tout;
      partialTrace(Tout, T, la, lb, on);
      T = Tout;

    }

  template<typename uni10_type>
    void partialTrace(UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& Tin, uni10_int32 la, uni10_int32 lb, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg(!((*Tin.status) & Tin.HAVEELEM), "%s" ,"Cannot perform contraction of two tensors before setting their elements.");
      uni10_error_msg(!(Tin.bonds->size() > 2), "%s", "The number of bonds must larger than 2 for performing partialTrace.");
      uni10_error_msg(&Tout == &Tin, "%s", 
          "The address of Tin and Tout are the same. Please use void partialTrace(UniTensor<uni10_type>& T, uni10_int32 la, uni10_int32 lb, UNI10_INPLACE on) instead.");


      uni10_int32 bondNum = Tin.bonds->size();
      std::vector<Bond> newBonds;
      std::vector<uni10_int32>newLabels(bondNum - 2, 0);
      std::vector<uni10_int32>rsp_labels(bondNum);
      uni10_int32 ia, ib;
      uni10_int32 enc = 0;

      for(uni10_uint64 l = 0; l < Tin.labels->size(); l++){
        if((*Tin.labels)[l] == la)
          ia = l;
        else if((*Tin.labels)[l] == lb)
          ib = l;
        else{
          newBonds.push_back((*Tin.bonds)[l]);
          newLabels[enc] = (*Tin.labels)[l];
          rsp_labels[enc] = (*Tin.labels)[l];
          enc++;
        }
      }

      uni10_error_msg(!(enc == newLabels.size()), "%s", "Cannot find the two bonds with the given two labels.");

      Tout.assign(newBonds);
      Tout.setLabel(newLabels);
      rsp_labels[bondNum - 2] = (*Tin.labels)[ia];
      rsp_labels[bondNum - 1] = (*Tin.labels)[ib];
      ia = bondNum - 2;
      ib = bondNum - 1;

      UniTensor<uni10_type> pTin;
      permute(pTin, Tin, rsp_labels, *Tout.RBondNum, INPLACE);

      tensor_tools::traceByRow(Tout.paras, pTin.paras, ia, ib, pTin.style);

      *Tout.status |= UniTensor<uni10_type>::GET_HAVEELEM();

    }

};

#endif
