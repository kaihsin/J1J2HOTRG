#ifndef __UNI10_LINALG_INPLACE_QL_H__
#define __UNI10_LINALG_INPLACE_QL_H__

#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  template<typename uni10_type>
    void ql( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& Q, Matrix<uni10_type>& L, UNI10_INPLACE on  );

  template<typename uni10_type>
    void ql( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& Q, Matrix<uni10_type>& L, UNI10_INPLACE on  ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum < Mij.Cnum, "%s", "Cannot perform QL decomposition when Rnum < Cnum. Nothing to do." );

      Q.assign(Mij.Rnum, Mij.Cnum);
      L.assign(Mij.Cnum, Mij.Cnum);

      matrixQL(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &Q.elem, &L.elem);

    }

}

#endif
