#ifndef __UNI10_LINALG_INPLACE_EIGH_H__
#define __UNI10_LINALG_INPLACE_EIGH_H__

#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  template<typename uni10_type>
    void eigh( const Matrix<uni10_type>& Mij, Matrix<uni10_double64>& Eig, Matrix<uni10_type>& EigVec, UNI10_INPLACE on );

  template<typename uni10_type>
    void eigh( const Matrix<uni10_type>& Mij, Matrix<uni10_double64>& Eig, Matrix<uni10_type>& EigVec, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg(Mij.Rnum != Mij.Cnum, "%s", "Setting a wrong flag of uni10_Inplace." );

      //GPU_NOT_READY
      Eig.assign(Mij.Rnum, Mij.Cnum, true);
      EigVec.assign(Mij.Rnum, Mij.Cnum);

      matrixEigh(&Mij.elem, &Mij.diag, &Mij.Cnum, &Eig.elem, &EigVec.elem);

    }

}

#endif
