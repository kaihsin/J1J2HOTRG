#ifndef __UNI10_LINALG_INPLACE_OTIMESARGS_H__
#define __UNI10_LINALG_INPLACE_OTIMESARGS_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_otimess.h"

namespace uni10{

  // DOES NOT OPTIMIZE.
  template<typename Res> 
    void _otimes_args(Res& TMout, std::vector<std::pair<void*, int> >& umlist);

  template<typename Res, typename Obj, typename... Args> 
    void _otimes_args(Res& TMout, std::vector<std::pair<void*, int> >& umlist, const Obj& _m2, const Args&... args);

  template<typename Res, typename Obj, typename... Args> 
    void otimes_args(Res& TMout, const Obj& _m1, const Args&... args);


  template<typename Res> 
    void _otimes_args(Res& TMout, std::vector<std::pair<void*, int> >& umlist){

      bool is_pure = true;
      for(int i = 0; i < umlist.size()-1; i++)
        if(umlist[i].second != umlist[i+1].second){
          is_pure = false;
          break;
        }

      if(is_pure)
        otimess_pure(TMout, umlist, INPLACE);
      else{
        otimess_mix(TMout, umlist, INPLACE);
      }

    }

  template<typename Res, typename Obj, typename... Args> 
    void _otimes_args(Res& TMout, std::vector<std::pair<void*, int> >& umlist, const Obj& _m1, const Args&... args) {

      umlist.push_back(std::pair<void*, int>((void*)&_m1, _m1.typeID()));
      _otimes_args(TMout, umlist, args...);

    }

  template<typename Res, typename Obj, typename... Args> 
    void otimes_args(Res& TMout, const Obj& _m1, const Args&... args) {

      std::vector<std::pair<void*, int> > umlist;
      umlist.push_back(std::pair<void*, int>((void*)&_m1, _m1.typeID()));
      _otimes_args(TMout, umlist, args...);

    }

}

#endif
