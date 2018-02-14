#ifndef __UNI10_LINALG_RQ_H__
#define __UNI10_LINALG_RQ_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"
//#include "uni10/uni10_api/linalg_typemix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > rq( const Block<uni10_type>& Mij );

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > rq( const Block<uni10_type>& Mij ){

      uni10_error_msg(Mij.Rnum > Mij.Cnum, "%s", "Cannot perform RQ decomposition when Rnum > Cnum. Nothing to do." );

      std::vector<Matrix<uni10_type> > outs;
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Rnum));
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, Mij.Cnum));

      matrixRQ(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem);

      return outs;

    }

}

#endif
