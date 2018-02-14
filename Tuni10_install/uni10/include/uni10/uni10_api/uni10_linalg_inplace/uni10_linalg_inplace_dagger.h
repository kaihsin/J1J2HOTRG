#ifndef __UNI10_LINALG_INPLACE_DAGGER_H__
#define __UNI10_LINALG_INPLACE_DAGGER_H__

#include "uni10/uni10_api/Matrix.h"

namespace uni10{


  template<typename uni10_type>
    void dagger( Matrix<uni10_type>& Mij, const Matrix<uni10_type>& ori_Mij, UNI10_INPLACE on );

  template<typename uni10_type>
    void dagger( Matrix<uni10_type>& Mij, UNI10_INPLACE on );

  template<typename uni10_type>
    void dagger( Matrix<uni10_type>& Mij, const Matrix<uni10_type>& ori_Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      Mij.assign(ori_Mij.Cnum, ori_Mij.Rnum, ori_Mij.diag);
      setDagger(&ori_Mij.elem, &ori_Mij.Rnum, &ori_Mij.Cnum, &Mij.elem);

    }

  template<typename uni10_type>
    void dagger( Matrix<uni10_type>& Mij, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      setDagger(&Mij.elem, &Mij.Rnum, &Mij.Cnum);

      uni10_uint64 tmp = Mij.Rnum;
      Mij.Rnum = Mij.Cnum;
      Mij.Cnum = tmp;

    }

}

#endif
