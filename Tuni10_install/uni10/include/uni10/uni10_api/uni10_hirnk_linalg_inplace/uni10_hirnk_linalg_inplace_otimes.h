#ifndef __UNI10_HIGH_RANK_LINALG_INPLACE_OTIMES_H__
#define __UNI10_HIGH_RANK_LINALG_INPLACE_OTIMES_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_contract.h"

namespace uni10{

  template<typename To, typename T, typename U>
    void otimes( UniTensor<To>& Tc, const UniTensor<T>& Ta, const UniTensor<U>& Tb, UNI10_INPLACE on);

  template<typename To, typename T, typename U>
    void otimes( Matrix<To>& Mc, const Block<T>& Ta, const Block<U>& Tb, UNI10_INPLACE on);

  template<typename To, typename T, typename U>
    void otimes( UniTensor<To>& Tc, const UniTensor<T>& Ta, const UniTensor<U>& Tb, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      
      UniTensor<T> T1 = Ta;
      UniTensor<U> T2 = Tb;
      std::vector<uni10_int> label1(T1.bondNum());
      std::vector<uni10_int> label2(T2.bondNum());
      for(uni10_uint64 i = 0; i < T1.bondNum(); i++){
        if(i < T1.inBondNum())
          label1[i] = i;
        else
          label1[i] = T2.inBondNum() + i;
      }
      for(uni10_uint64 i = 0; i < T2.bondNum(); i++){
        if(i < T2.inBondNum())
          label2[i] = i + T1.inBondNum();
        else
          label2[i] = i + T1.bondNum();
      }
      T1.setLabel(label1);
      T2.setLabel(label2);

      contract(Tc, T1, T2, true, on);

    }

  template<typename To, typename T, typename U>
    void otimes( Matrix<To>& Mc, const Block<T>& Ta, const Block<U>& Tb, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      UniTensor<T> T1(Ta);
      UniTensor<U> T2(Tb);
      std::vector<uni10_int> label1(T1.bondNum());
      std::vector<uni10_int> label2(T2.bondNum());
      for(uni10_uint64 i = 0; i < T1.bondNum(); i++){
        if(i < T1.inBondNum())
          label1[i] = i;
        else
          label1[i] = T2.inBondNum() + i;
      }
      for(uni10_uint64 i = 0; i < T2.bondNum(); i++){
        if(i < T2.inBondNum())
          label2[i] = i + T1.inBondNum();
        else
          label2[i] = i + T1.bondNum();
      }
      T1.setLabel(label1);
      T2.setLabel(label2);

      UniTensor<To> Tc;
      contract(Tc, T1, T2, true, INPLACE);
      Mc = Tc.const_getBlock();

    }

};

#endif
