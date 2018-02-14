#ifndef __UNI10_LINALG_GETDIAG_H__
#define __UNI10_LINALG_GETDIAG_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"
//#include "uni10/uni10_api/linalg_typemix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type> 
    Matrix<uni10_type> getDiag( const Block<uni10_type>& A );

  template<typename uni10_type> 
    Matrix<uni10_type> getDiag( const Block<uni10_type>& A ){
      if(A.diag){
        Matrix<uni10_type> D(A);
        return D;
      }
      else{
        Matrix<uni10_type> D(A.Rnum, A.Cnum, true);
        uni10_error_msg(true, "%s", "Developping!!!\n");
        return D;
      }
    }

}

#endif
