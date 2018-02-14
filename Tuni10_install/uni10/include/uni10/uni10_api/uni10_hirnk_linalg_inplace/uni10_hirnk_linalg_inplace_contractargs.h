#ifndef __UNI10_HIGH_RANK_LINALG_INPLACE_CONTRACT_ARGS_H__
#define __UNI10_HIGH_RANK_LINALG_INPLACE_CONTRACT_ARGS_H__

#include "uni10/uni10_api/UniTensor.h"
#include "uni10/uni10_api/Network.h"
#include "uni10/uni10_api/Network_dev.h"

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_contracts.h"

namespace uni10{

  //
  // contract_args with Network_dev.
  // 
  template<typename UniTo>
    void _contract_args( UniTo& Tout, Network_dev& ctr_net, uni10_uint64& idx);

  template<typename UniTo, typename UniTi, typename ... Args>
    void _contract_args( UniTo& Tout, Network_dev& ctr_net, uni10_uint64& idx, const UniTi& T1, const Args&... args);

  template<typename UniTo, typename UniTi, typename ... Args>
    void contract_args( UniTo& Tout, Network_dev& ctr_net, const UniTi& T1, const Args&... args);


  template<typename UniTo>
    void _contract_args( UniTo& Tout, Network_dev& ctr_net, uni10_uint64& idx){
      ctr_net.launch(Tout);
    }

  template<typename UniTo, typename UniTi, typename ... Args>
    void _contract_args( UniTo& Tout, Network_dev& ctr_net, uni10_uint64& idx, const UniTi& T1, const Args&... args){
      ctr_net.putTensor(idx, T1);
      idx++;
      _contract_args(Tout, ctr_net, idx, args...);
    }

  template<typename UniTo, typename UniTi, typename ... Args>
    void contract_args( UniTo& Tout, Network_dev& ctr_net, const UniTi& T1, const Args&... args){
      uni10_uint64 idx = 0;
      _contract_args(Tout, ctr_net, idx, T1, args...);
    }

  //
  // contract_args without Network. Before using this function, users have to set the suitable labers on each UniTensors' bonds.
  // 
  template<typename UniTo>
    void _contract_args( UniTo& Tout, bool fast, std::vector< std::pair<void*, int> >& ulist);

  template<typename UniTo, typename UniTi, typename ... Args>
    void _contract_args( UniTo& Tout, bool fast, std::vector< std::pair<void*, int> >& ulist, const UniTi& T1, const Args&... args);

  template<typename UniTo, typename UniTi, typename ... Args>
    void contract_args( UniTo& Tout, bool fast, const UniTi& T1, const Args&... args);

  template<typename UniTo>
    void _contract_args( UniTo& Tout, bool fast, std::vector< std::pair<void*, int> >& ulist){

      std::pair<bool, uni10_int> isMix_TrLen(false, 0);

      bool isTc = false;

      uni10_int cnt_Tr = 0;

      for(uni10_uint64 i = 0; i < ulist.size(); i++){

        if(ulist[i].second == 2)
          isTc = true;

        if(ulist[i].second == 1 && !isTc)
          cnt_Tr++;
       
        if(i == 0)
          continue;

        else if(ulist[i-1].second != ulist[i].second){

          isMix_TrLen = std::pair<bool, uni10_int>(true, cnt_Tr);

        }

      }

      if(isMix_TrLen.first == false){

        contracts_pure(Tout, ulist, fast, INPLACE);

      }
      else{

        contracts_mix(Tout, ulist, isMix_TrLen, fast, INPLACE);
      }

    }

  template<typename UniTo, typename UniTi, typename ... Args>
    void _contract_args( UniTo& Tout, bool fast, std::vector< std::pair<void*, int> >& ulist, const UniTi& T1, const Args&... args){

      ulist.push_back(std::pair<void*, int>((void*)&T1, T1.typeID()));
      _contract_args(Tout, fast, ulist, args...);

    }


  template<typename UniTo, typename UniTi, typename ... Args>
    void contract_args( UniTo& Tout, bool fast, const UniTi& T1, const Args&... args){

      std::vector< std::pair<void*, int> > ulist;
      _contract_args(Tout, fast, ulist, T1, args...);

    }

  //
  // contract_args with previous Network class, Network<T>.
  // After doing more restrictively test to the Network_dev, we will discard these functions.
  // 
  template<typename UniTo, typename U>
    void contract_args( UniTo& Tout, Network<U>& Ta, uni10_uint64& idx);

  template<typename UniTo, typename UniTi, typename U, typename ... Args>
    void contract_args( UniTo& Tout, Network<U>& ctr_net, const UniTi& T1, const Args&... args);

  template<typename UniTo, typename UniTi, typename U, typename ... Args>
    void contract_args( UniTo& Tout, Network<U>& ctr_net, uni10_uint64& idx, const UniTi& T1, const Args&... args);


  template<typename UniTo, typename U>
    void contract_args( UniTo& Tout, Network<U>& ctr_net, uni10_uint64& idx){
      ctr_net.launch(Tout);
    }

  template<typename UniTo, typename UniTi, typename U, typename ... Args>
    void contract_args( UniTo& Tout, Network<U>& ctr_net, const UniTi& T1, const Args&... args){
      uni10_uint64 idx = 0;
      contract_args(Tout, ctr_net, idx, T1, args...);
    }

  template<typename UniTo, typename UniTi, typename U, typename ... Args>
    void contract_args( UniTo& Tout, Network<U>& ctr_net, uni10_uint64& idx, const UniTi& T1, const Args&... args){
      ctr_net.putTensor(idx, T1);
      idx++;
      contract_args(Tout, ctr_net, idx, args...);
    }

};

#endif
