#ifndef __UNI10_LINALG_CONJ_H__
#define __UNI10_LINALG_CONJ_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"
//#include "uni10/uni10_api/linalg_typemix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type>
    Matrix<uni10_type> conj( const Block<uni10_type>& Mij );

  template<typename uni10_type>
    Matrix<uni10_type> conj( const Block<uni10_type>& Mij ){

      Matrix<uni10_type> MijConj(Mij.Rnum, Mij.Cnum, Mij.diag);
      setConjugate(&Mij.elem, &Mij.elem.__elemNum, &MijConj.elem);

      return MijConj;

    }

}

#endif
