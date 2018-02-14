#ifndef __UNI10_LINALG_INPLACE_RESIZE_H__
#define __UNI10_LINALG_INPLACE_RESIZE_H__

#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  template<typename _uni10_type> 
    void resize( Matrix<_uni10_type>& M , uni10_uint64 row, uni10_uint64 col, UNI10_INPLACE on);

  template<typename _uni10_type> 
    void resize( Matrix<_uni10_type>& Mout , const Matrix<_uni10_type>& Min, uni10_uint64 row, uni10_uint64 col, UNI10_INPLACE on);

  template<typename _uni10_type> 
    void resize( Matrix<_uni10_type>& M , uni10_uint64 row, uni10_uint64 col, UNI10_INPLACE on){ 

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      Matrix<_uni10_type> Mout;
      resize(Mout, M, row,  col, on);
      M = Mout;

    }

  template<typename _uni10_type> 
    void resize( Matrix<_uni10_type>& Mout , const Matrix<_uni10_type>& Min, uni10_uint64 row, uni10_uint64 col, UNI10_INPLACE on){ 

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      Mout.assign(row, col, Min.diag);
      uni10_bool fixHead = true;
      resize(Mout.elem, Mout.Rnum, Mout.Cnum, Min.elem, Min.Rnum, Min.Cnum, fixHead);

    }

}

#endif
