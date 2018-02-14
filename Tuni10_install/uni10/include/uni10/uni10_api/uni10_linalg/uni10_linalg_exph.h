#ifndef __UNI10_LINALG_EXPH_H__
#define __UNI10_LINALG_EXPH_H__

#include "uni10/uni10_api/uni10_linalg_inplace/uni10_linalg_inplace_dot.h"

namespace uni10{

  template<typename uni10_type>
    Matrix<uni10_type> exph( uni10_double64 a, const Block<uni10_type>& mat);

  template<typename uni10_type>
    Matrix<uni10_type> exph( uni10_double64 a, const Block<uni10_type>& mat){

      std::vector< Matrix<uni10_type> > rets = eigh( mat );

      Matrix<uni10_type> UT, EXPT;
      dagger(UT, rets[1], INPLACE);

      vectorExp( &a, &rets[0].elem, &rets[0].Rnum );

      dot_args(EXPT, UT, rets[0], rets[1]);

      return EXPT;

    }

}

#endif
