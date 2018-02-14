#ifndef __UNI10_LINALG_INPLACE_SVD_H__
#define __UNI10_LINALG_INPLACE_SVD_H__

#include "uni10/uni10_api/Matrix.h"

namespace uni10{


  template<typename uni10_type>
    void svd( const Block<uni10_type>& Mij, Matrix<uni10_type>& U, Matrix<uni10_type>& S, Matrix<uni10_type>& vT, UNI10_INPLACE on );

  template<typename uni10_type>
    void svd( const Block<uni10_type>& Mij, Matrix<uni10_type>& U, Matrix<uni10_type>& S, Matrix<uni10_type>& vT, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_uint64 min = Mij.Rnum < Mij.Cnum ? Mij.Rnum : Mij.Cnum;      //min = min(Rnum,Cnum)
      S.assign(min, min, true);

      UELEM(uni10_elem, _package, _type)<uni10_type>* U_elem  = NULL;     // pointer to a real matrix
      UELEM(uni10_elem, _package, _type)<uni10_type>* vT_elem = NULL;     // pointer to a real matrix

      if(&U != NULL){
        U.assign(Mij.Rnum, min);
        U_elem = &U.elem;
      }

      if(&vT != NULL){
        vT.assign(min, Mij.Cnum);
        vT_elem = &vT.elem;
      }

      matrixSVD(&Mij.elem, &Mij.diag, &Mij.Rnum, &Mij.Cnum, U_elem, &S.elem, vT_elem);

    }

}

#endif
