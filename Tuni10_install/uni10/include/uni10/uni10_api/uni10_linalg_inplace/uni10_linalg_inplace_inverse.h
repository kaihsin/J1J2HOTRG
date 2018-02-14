#ifndef __UNI10_LINALG_INPLACE_INVERSE_H__
#define __UNI10_LINALG_INPLACE_INVERSE_H__

#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  template<typename uni10_type>
    void inverse( Matrix<uni10_type>& Mij, UNI10_INPLACE on );

  template<typename uni10_type>
    void inverse( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(!(Mij.Rnum == Mij.Cnum), "%s", "Cannot perform inversion on a non-square matrix." );

      matrixInv(&Mij.elem, &Mij.Rnum, &Mij.diag);

    }

}

#endif
