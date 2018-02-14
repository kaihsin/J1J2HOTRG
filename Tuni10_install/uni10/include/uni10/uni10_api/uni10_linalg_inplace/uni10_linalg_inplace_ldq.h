#ifndef __UNI10_LINALG_INPLACE_LDQ_H__
#define __UNI10_LINALG_INPLACE_LDQ_H__

#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  template<typename uni10_type>
    void ldq( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& L, Matrix<uni10_type>& D, Matrix<uni10_type>& Q, UNI10_INPLACE on  );

  template<typename uni10_type>
    void ldq( const Matrix<uni10_type>& Mij, Matrix<uni10_type>& L, Matrix<uni10_type>& D, Matrix<uni10_type>& Q, UNI10_INPLACE on  ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "%s", "Cannot perform LDQ decomposition when Rnum > Cnum. Nothing to do." );

      L.assign(Mij.Rnum, Mij.Rnum);
      D.assign(Mij.Rnum, Mij.Rnum, true);
      Q.assign(Mij.Rnum, Mij.Cnum);

      matrixLDQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &L.elem, &D.elem, &Q.elem);

    }

}

#endif
