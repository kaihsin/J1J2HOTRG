#ifndef __UNI10_LINALG_INPLACE_DOTARGS_H__
#define __UNI10_LINALG_INPLACE_DOTARGS_H__

#include "uni10/uni10_api/Matrix.h"
#include "uni10_linalg_inplace_dots.h"

namespace uni10{

  // DOES NOT OPTIMIZE.
  template<typename Res> 
    void _dot_args(Res& mout, std::vector<std::pair<void*, int> >& mlist);

  template<typename Res, typename Obj, typename... Args> 
    void _dot_args(Res& mout, std::vector<std::pair<void*, int> >& mlist, const Obj& _m2, const Args&... args);

  template<typename Res, typename Obj, typename... Args> 
    void dot_args(Res& mout, const Obj& _m1, const Args&... args);



  template<typename Res> 
    void _dot_args(Res& mout, std::vector<std::pair<void*, int> >& mlist){

      bool is_pure = true;
      for(int i = 0; i < mlist.size()-1; i++)
        if(mlist[i].second != mlist[i+1].second){
          is_pure = false;
          break;
        }

      if(is_pure)
        dots_pure(mout, mlist, INPLACE);
      else{
        dots_mix(mout, mlist, INPLACE);
      }

    }

  template<typename Res, typename Obj, typename... Args> 
    void _dot_args(Res& mout, std::vector<std::pair<void*, int> >& mlist, const Obj& _m1, const Args&... args) {

      mlist.push_back(std::pair<void*, int>((void*)&_m1, _m1.typeID()));
      _dot_args(mout, mlist, args...);

    }

  template<typename Res, typename Obj, typename... Args> 
    void dot_args(Res& mout, const Obj& _m1, const Args&... args) {

      std::vector<std::pair<void*, int> > mlist;
      mlist.push_back(std::pair<void*, int>((void*)&_m1, _m1.typeID()));
      _dot_args(mout, mlist, args...);
      /*
      if(_m1.elemNum() == 0)
        _m1 = _m2;
      else
        dot(_m1, _m2, INPLACE);
      dot_args(_m1, args...);
      */
    }

}

#endif
