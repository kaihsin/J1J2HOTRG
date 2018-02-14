#ifndef __UNI10_HIGH_RANK_LINALG_INPLACE_CONTRACTS_H__
#define __UNI10_HIGH_RANK_LINALG_INPLACE_CONTRACTS_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_contract.h"

namespace uni10{

  // The driver fucntions of Non-type UniTensor contraction.
  inline void contract_dd(void* m3, void* m1, void* m2, bool fast){

    contract( *((UniTensor<uni10_double64>*)m3), *((UniTensor<uni10_double64>*)m1), *((UniTensor<uni10_double64>*)m2), fast, INPLACE);

  }

  inline void contract_dz(void* m3, void* m1, void* m2, bool fast){

    contract( *((UniTensor<uni10_complex128 >*)m3), *((UniTensor<uni10_double64>*)m1), *((UniTensor<uni10_complex128 >*)m2), fast, INPLACE);

  }

  inline void contract_zd(void* m3, void* m1, void* m2, bool fast){

    contract( *((UniTensor<uni10_complex128 >*)m3), *((UniTensor<uni10_complex128 >*)m1), *((UniTensor<uni10_double64>*)m2), fast, INPLACE);

  }

  inline void contract_zz(void* m3, void* m1, void* m2, bool fast){

    contract( *((UniTensor<uni10_complex128 >*)m3), *((UniTensor<uni10_complex128 >*)m1), *((UniTensor<uni10_complex128 >*)m2), fast, INPLACE);

  }

  // Function pointers of vector of contracts' driver functions.
  static void (*contract_driver[])(void* m3, void* m1, void* m2, bool fast) = {contract_dd, contract_zd, contract_dz, contract_zz};


  // contract with mix type.
  template<typename T>
    void contracts_mix(UniTensor<T>& mout, std::vector< std::pair<void*, int> >& _ulist, bool fast, UNI10_INPLACE on);

  template<typename T> 
    void contracts_pure(UniTensor<T>& mout, std::vector< std::pair<void*, int> >& _ulist, bool fast, UNI10_INPLACE on);

  // Uni10_API level.
  template<typename T> 
    void contracts(UniTensor<T>& mout, std::vector< UniTensor<T>* >& _ulist, std::pair<bool, uni10_int>& isMix_TrLen, bool fast, UNI10_INPLACE on);

  // contract with mix type.
  template<typename T>
    void contracts_mix(UniTensor<T>& _Tout, std::vector< std::pair<void*, int> >& _ulist, std::pair<bool, uni10_int>& isMix_TrLen, bool fast, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace");
      uni10_error_msg(_ulist.size() < 2, "%s", "Too few UniTensors in the contraction list [N >= 2]. ");
      uni10_error_msg(_Tout.typeID() == 1, "%s", 
          "There are complex matrices in the input arguments. Hence, the type of output matrix must be complex.");

      //
      // Check Tout is in the ulist or not.
      //
      
      bool inlist = false;
      for(uni10_uint64 i = 0; i < _ulist.size(); i++)
        if(_ulist[i].first == &_Tout){
          inlist = true;
          break;
        }

      UniTensor<T>* Tout = NULL;

      if(inlist){
        Tout = new UniTensor<T>(UniTensor<T>());
      }else{
        Tout = &_Tout;
      }

      uni10_uint64 complex1st = isMix_TrLen.second;

      if(isMix_TrLen.second > 1){

        UniTensor<uni10_double64> Tr1, Tr2;

        contract_driver[0]((void*)&Tr1, _ulist[0].first, _ulist[1].first, fast);

        uni10_int i = 2;

        for( i = 2; i < isMix_TrLen.second; i++){

          Tr2 = *(UniTensor<uni10_double64>*)_ulist[i].first;

          if(i % 2 == 0)
            contract_driver[0]((void*)&Tr2, (void*)&Tr1, _ulist[i].first, fast);
          else
            contract_driver[0]((void*)&Tr1, (void*)&Tr2, _ulist[i].first, fast);

        }

        if(i % 2 == 0){
          contract_driver[2]((void*)Tout, (void*)&Tr1, _ulist[complex1st].first, fast);
        }

        if(i % 2 == 1){
          contract_driver[2]((void*)Tout, (void*)&Tr2, _ulist[complex1st].first, fast);
        }

      }else{

        uni10_int rest_idx = isMix_TrLen.second == 0 ? 1 : 0;

        if(_ulist[rest_idx].second == 1){

          contract_driver[1]((void*)Tout, _ulist[complex1st].first, _ulist[rest_idx].first, fast);

        }else{

          contract_driver[3]((void*)Tout, _ulist[complex1st].first, _ulist[rest_idx].first, fast);

        }

      }

      uni10_uint64 offset = (isMix_TrLen.second > 1) ? isMix_TrLen.second + 1 : 2;
      uni10_uint64 j = 0;

      UniTensor<uni10_complex128> tmpT1;

      for( ; j < _ulist.size()-offset; j++){

        int driver_type = _ulist[j+offset].second == 1 ? 1 : 3;

        if(j % 2 == 0)
          contract_driver[driver_type]((void*)&tmpT1, (void*)Tout, _ulist[j+offset].first, fast);
        else
          contract_driver[driver_type]((void*)Tout, (void*)&tmpT1, _ulist[j+offset].first, fast);

      }

      if(j % 2 == 1){

        *Tout = tmpT1;

      }

      if(inlist){
        _Tout = *Tout;
        delete Tout;
      }

    }


  template<typename T> 
    void contracts_pure(UniTensor<T>& _Tout, std::vector< std::pair<void*, int> >& _ulist, bool fast, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace");
      uni10_error_msg(_Tout.typeID() != _ulist[0].second, "%s", "Set a wrong type to the output UniTensor.");
      
      //
      // Check Tout is in the ulist or not.
      //
      bool inlist = false;

      for(uni10_uint64 i = 0; i < _ulist.size(); i++)
        if(_ulist[i].first == &_Tout){
          inlist = true;
          break;
        }

      UniTensor<T>* Tout = NULL;

      if(inlist){
        Tout = new UniTensor<T>(UniTensor<T>());
      }else{
        Tout = &_Tout;
      }
           
      int driver_type = (Tout->typeID() == 1) ? 0 : 3;

      contract_driver[driver_type]((void*)Tout, _ulist[0].first, _ulist[1].first, fast);

      UniTensor<T> tmpT1;

      uni10_uint64 i = 2;

      for( ; i < _ulist.size(); i++){

        if(i % 2 == 0)
          contract_driver[driver_type]((void*)&tmpT1, (void*)Tout, _ulist[i].first, fast);
        else
          contract_driver[driver_type]((void*)Tout, (void*)&tmpT1, _ulist[i].first, fast);

      }

      if(i % 2 == 1){

        *Tout = tmpT1;

      }

      if(inlist){
        _Tout = *Tout;
        delete Tout;
      }

    };

  template<typename T> 
    void contracts(UniTensor<T>& Tout, std::vector< UniTensor<T>* >& _ulist, bool fast, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace");

      contract(Tout, *_ulist[0], *_ulist[1], fast, on);

      uni10_uint64 i;

      UniTensor<T> tmpT1;

      for( i = 2; i < _ulist.size(); i++){

        if(i % 2 == 0)
          contract(tmpT1, Tout, *_ulist[i], fast, on);
        else
          contract(Tout, tmpT1, *_ulist[i], fast, on);

      }

      if(i % 2 == 1){

        Tout = tmpT1;

      }

    }

}; /* End of namespace. */

#endif
