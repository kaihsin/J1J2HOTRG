#ifndef __UNI10_LINALG_NORM_H__
#define __UNI10_LINALG_NORM_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"
//#include "uni10/uni10_api/linalg_typemix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type>
    uni10_double64 norm( const Block<uni10_type>& Mij );

  template<typename uni10_type>
    uni10_double64 norm( const Block<uni10_type>& Mij ){

      uni10_int32 inc = 1;
      uni10_uint64 elemNum = Mij.diag ? std::min(Mij.Rnum, Mij.Cnum) : Mij.Rnum * Mij.Cnum;
      return vectorNorm(&Mij.elem, &elemNum, &inc);

    }

}

#endif
