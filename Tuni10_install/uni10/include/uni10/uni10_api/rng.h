#ifndef __UNI10_RNG_H__
#define __UNI10_RNG_H__

#include <chrono>

#include "uni10/uni10_elem_rng.h"
#include "uni10/uni10_api/linalg_inplace.h"

#define uni10_rand(M, eng, dis, up, dn, seed)\
  do{ \
    uni10_int32  seed1 = seed;\
    uni10_uint64 num = M.elemNum();\
    if(seed1 == -1)\
      seed1 = std::chrono::system_clock::now().time_since_epoch().count();\
    uni10_elem_rng(M, num, eng, dis, up, dn, seed1);\
  }while(0);

#define uni10_orthoRand(M, eng, dis, up, dn, seed)\
  do{ \
    uni10_int32 seed1 = seed;\
    uni10_uint64 num = M.elemNum();\
    if(seed1 == -1)\
      seed1 = std::chrono::system_clock::now().time_since_epoch().count();\
    uni10_elem_rng(M, num, eng, dis, up, dn, seed1);\
    M = uni10::svd(M)[0];\
  }while(0);


/*)
namespace uni10{

  template <typename uni10_type>
    void identity(Matrix<uni10_type>& M){
      
      uni10_error_msg(M.elemNum() == 0, "%s", "The matrix has not been initialized!!!");
      M.diag = true;
      if(! M.elem.empty())
        M.elem.clear();
      M.elem.init(M.Rnum, M.Cnum, M.diag);
      for(uni10_uint64 i = 0; i < M.elem.__elemNum; i++)
        M.elem.__elem[i] = 1.;
    }


}
*/

#endif
