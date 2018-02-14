#ifndef __UNI10_LINALG_DOT_H__
#define __UNI10_LINALG_DOT_H__

#include "uni10/uni10_api/uni10_linalg_inplace/uni10_linalg_inplace_dot.h"

namespace uni10{

  template<typename T> 
    Matrix<uni10_complex128> dot( const Block<uni10_complex128>& A, const Block<T>& B );

  template<typename T> 
    Matrix<T> dot( const Block<uni10_double64>& A, const Block<T>& B );

  template<typename T> 
    Matrix<uni10_complex128> dot( const Block<uni10_complex128>& A, const Block<T>& B ){

      Matrix<uni10_complex128> C;
      dot(C, A, B, INPLACE);
      return C;

    }

  template<typename T> 
    Matrix<T> dot( const Block<uni10_double64>& A, const Block<T>& B ){

      Matrix<T> C;
      dot(C, A, B, INPLACE);
      return C;

    }

}

#endif
