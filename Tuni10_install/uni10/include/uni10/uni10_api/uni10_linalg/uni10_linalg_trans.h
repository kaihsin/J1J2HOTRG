#ifndef __UNI10_LINALG_TRANS_H__
#define __UNI10_LINALG_TRANS_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"
//#include "uni10/uni10_api/linalg_typemix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type>
    Matrix<uni10_type> transpose( const Block<uni10_type>& Mij ){

      Matrix<uni10_type> MijT(Mij.Cnum, Mij.Rnum, Mij.diag);
      setTranspose(&Mij.elem, &Mij.Rnum, &Mij.Cnum, &MijT.elem);

      return MijT;

    }

}

#endif
