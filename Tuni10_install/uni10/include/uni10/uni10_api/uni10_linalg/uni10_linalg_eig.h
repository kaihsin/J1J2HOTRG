#ifndef __UNI10_LINALG_EIG_H__
#define __UNI10_LINALG_EIG_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type>
    std::vector< Matrix<uni10_complex128> > eig( const Block<uni10_type>& Mij );

  template<typename uni10_type>
    std::vector< Matrix<uni10_complex128> > eig( const Block<uni10_type>& Mij ){

      std::vector< Matrix<uni10_complex128> > outs;
      outs.push_back(Matrix<uni10_complex128>(Mij.Rnum, Mij.Cnum, true));
      outs.push_back(Matrix<uni10_complex128>(Mij.Rnum, Mij.Cnum));
      matrixEig(&Mij.elem, &Mij.diag, &Mij.Cnum, &outs[0].elem, &outs[1].elem);

      return outs;

    }

}

#endif
