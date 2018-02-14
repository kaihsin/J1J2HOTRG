#ifndef __UNI10_LINALG_INPLACE_DOT_H__
#define __UNI10_LINALG_INPLACE_DOT_H__

#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  void dot( Matrix<uni10_double64>& A, const Block<uni10_double64>& B, UNI10_INPLACE on );

  template<typename T> 
    void dot( Matrix<uni10_complex128>& A, const Block<T>& B, UNI10_INPLACE on );

  template<typename To, typename T, typename U> 
    void dot( Matrix<To>& C, const Block<T>& A, const Block<U>& B, UNI10_INPLACE on );

  template<typename T> 
    void dot( Matrix<uni10_complex128>& A, const Block<T>& B, UNI10_INPLACE on ){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(A.Cnum != B.Rnum, "%s", "The dimensions of the two matrices do not match for matrix multiplication.");

      Matrix<uni10_complex128> C(A.Rnum, B.Cnum, A.diag && B.diag);

      matrixDot(&A.elem, &A.diag, &B.elem, &B.diag, &A.Rnum, &B.Cnum, &A.Cnum, &C.elem);

      A = C;

    }

  template<typename To, typename T, typename U> 
    void dot( Matrix<To>& C, const Block<T>& A, const Block<U>& B, UNI10_INPLACE on ){
      
      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      uni10_error_msg(A.Cnum != B.Rnum, "%s", "The dimensions of the two matrices do not match for matrix multiplication.");
      //Lack of error msgs.
      //
      C.assign(A.Rnum, B.Cnum, A.diag && B.diag);

      matrixDot(&A.elem, &A.diag, &B.elem, &B.diag, &A.Rnum, &B.Cnum, &A.Cnum, &C.elem);

    }

}

#endif
