#ifndef __UNI10_LINALG_INVERSE_H__
#define __UNI10_LINALG_INVERSE_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"
//#include "uni10/uni10_api/linalg_typemix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type>
    Matrix<uni10_type> inverse( const Block<uni10_type>& Mij );

  template<typename uni10_type>
    Matrix<uni10_type> inverse( const Block<uni10_type>& Mij ){

      Matrix<uni10_type> invM(Mij);

      uni10_error_msg(!(Mij.Rnum == Mij.Cnum), "%s", "Cannot perform inversion on a non-square matrix." );

      matrixInv(&invM.elem, &Mij.Rnum, &Mij.diag);

      return invM;

    }

}

#endif
