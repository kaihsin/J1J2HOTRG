#ifndef __UNI10_LINALG_INPLACE_RQ_H__
#define __UNI10_LINALG_INPLACE_RQ_H__

#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  template<typename uni10_type>
    void rq( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& R, Matrix<uni10_type>& Q, UNI10_INPLACE on  );

  template<typename uni10_type>
    void rq( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& R, Matrix<uni10_type>& Q, UNI10_INPLACE on  ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "%s", "Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do." );

      R.assign(Mij.Rnum, Mij.Rnum);
      Q.assign(Mij.Rnum, Mij.Cnum);

      matrixRQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &R.elem, &Q.elem);

    }

}

#endif
