#ifndef __UNI10_HIGH_RANK_LINALG_INPLACE_COJ_H__
#define __UNI10_HIGH_RANK_LINALG_INPLACE_COJ_H__

#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{

  template<typename uni10_type>
    void conj( UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& Tin, UNI10_INPLACE on);

  template<typename uni10_type>
    void conj( UniTensor<uni10_type>& T, UNI10_INPLACE on);

  template<typename uni10_type>
    void conj( UniTensor<uni10_type>& Tout, const UniTensor<uni10_type>& Tin, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg(!(*Tin.status & UniTensor<uni10_type>::GET_HAVEBOND()), 
          "%s", "There is no bond in the tensor(scalar) to perform transposition.");

      Tout.assign(*Tin.bonds);
      Tout.setName(*Tin.name);
      Tout.setLabel(*Tin.labels);

      if((*Tin.status) & UniTensor<uni10_type>::GET_HAVEELEM())
        tensor_tools::conj(Tout.paras, Tin.paras, Tin.style);
        

      *Tout.status |= UniTensor<uni10_type>::GET_HAVEELEM();

    }

  template<typename uni10_type>
    void conj( UniTensor<uni10_type>& T, UNI10_INPLACE on){

      UniTensor<uni10_type> Tout;
      conj(Tout, T, on);
      T = Tout;

    }

};

#endif
