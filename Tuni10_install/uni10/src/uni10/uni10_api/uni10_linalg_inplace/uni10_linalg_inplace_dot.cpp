#include "uni10/uni10_api/uni10_linalg_inplace/uni10_linalg_inplace_dot.h"

namespace uni10{

  void dot( Matrix<uni10_double64>& A, const Block<uni10_double64>& B, UNI10_INPLACE on ){

    uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

    uni10_error_msg(A.Cnum != B.Rnum, "%s", "The dimensions of the two matrices do not match for matrix multiplication.");

    Matrix<uni10_double64> C(A.Rnum, B.Cnum, A.diag && B.diag);

    matrixDot(&A.elem, &A.diag, &B.elem, &B.diag, &A.Rnum, &B.Cnum, &A.Cnum, &C.elem);

    A = C;

  }

}

