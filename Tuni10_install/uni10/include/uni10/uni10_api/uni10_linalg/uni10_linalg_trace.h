#ifndef __UNI10_LINALG_TRACE_H__
#define __UNI10_LINALG_TRACE_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"
//#include "uni10/uni10_api/linalg_typemix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type>
    uni10_type trace( const Block<uni10_type>& Mij );

  template<typename uni10_type>
    uni10_type trace( const Block<uni10_type>& Mij ){

      uni10_type res = matrixTrace(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum);

      return res;

    }

}

#endif
