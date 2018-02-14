#ifndef __UNI10_LINALG_INPLACE_CONJ_H__
#define __UNI10_LINALG_INPLACE_CONJ_H__

#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  template<typename uni10_type>
    void conj( Matrix<uni10_type>& Mij, const Matrix<uni10_type>& ori_Mij, UNI10_INPLACE on );

  template<typename uni10_type>
    void conj( Matrix<uni10_type>& Mij, UNI10_INPLACE on );

  template<typename uni10_type>
    void conj( Matrix<uni10_type>& Mij, const Matrix<uni10_type>& ori_Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      Mij.assign(ori_Mij.Rnum, ori_Mij.Cnum);
      setConjugate(&ori_Mij.elem, &Mij.Rnum, &Mij.Cnum, &Mij.elem);

    }

  template<typename uni10_type>
    void conj( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      setConjugate(&Mij.elem, &Mij.elem.__elemNum);

    }

}

#endif
