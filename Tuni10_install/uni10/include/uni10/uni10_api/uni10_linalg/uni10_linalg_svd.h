#ifndef __UNI10_LINALG_SVD_H__
#define __UNI10_LINALG_SVD_H__

#include <vector>

#include "uni10/uni10_api/Block.h"
#include "uni10/uni10_api/Matrix.h"
//#include "uni10/uni10_api/linalg_typemix.h"

#if defined(CPU) && defined(LAPACK)
#include "uni10/uni10_elem_linalg.h"
#endif

namespace uni10{

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > svd( const Block<uni10_type>& Mij );

  template<typename uni10_type>
    std::vector< Matrix<uni10_type> > svd( const Block<uni10_type>& Mij ){

      std::vector<Matrix<uni10_type> > outs;

      uni10_uint64 min = Mij.Rnum < Mij.Cnum ? Mij.Rnum : Mij.Cnum;      //min = min(Rnum,Cnum)
      //GPU_NOT_READY
      outs.push_back(Matrix<uni10_type>(Mij.Rnum, min));
      outs.push_back(Matrix<uni10_type>(min, min, true));
      outs.push_back(Matrix<uni10_type>(min, Mij.Cnum));

      matrixSVD(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, &outs[0].elem, &outs[1].elem, &outs[2].elem);

      return outs;

    }


}

#endif
