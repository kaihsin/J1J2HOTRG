#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  void matrixDot(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_double64* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_linalg::matrixDot(A->__elem, B->__elem, *M, *N, *K, C->__elem);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy(C->__elem, B->__elem, B->__elemNum*sizeof(uni10_double64));
      uni10_linalg::diagRowMul(C->__elem, A->__elem, std::min(*M, *K), *N);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 data_col = std::min(*K, *N);

      for(int r = 0; r < (int)*M; r++)
        uni10_elem_copy(&C->__elem[r*(*N)], &A->__elem[r*(*K)], data_col *sizeof(uni10_double64));

      uni10_linalg::diagColMul(C->__elem, B->__elem, *M, data_col);

    }
    else{

      uni10_uint64 min = std::min(A->__elemNum, B->__elemNum);
      uni10_elem_copy(C->__elem, A->__elem, min*sizeof(uni10_double64));

      uni10_linalg::vectorMul(C->__elem, B->__elem, min);

    }

  }

  void matrixDot(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_linalg::matrixDot(A->__elem, B->__elem, *M, *N, *K, C->__elem);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy(C->__elem, B->__elem, B->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::diagRowMul(C->__elem, A->__elem, std::min(*M, *K), *N);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 data_col = std::min(*K, *N);

      for(int r = 0; r < (int)*M; r++)
        uni10_elem_copy(&C->__elem[r*(*N)], &A->__elem[r*(*K)], data_col*sizeof(uni10_complex128));

      uni10_linalg::diagColMul(C->__elem, B->__elem, *M, data_col);

    }
    else{

      uni10_uint64 min = std::min(A->__elemNum, B->__elemNum);
      uni10_elem_copy(C->__elem, A->__elem, min*sizeof(uni10_complex128));

      uni10_linalg::vectorMul(C->__elem, B->__elem, min);

    }

  }

  void matrixDot(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_linalg::matrixDot(A->__elem, B->__elem, *M, *N, *K, C->__elem);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_copy(C->__elem, B->__elem, B->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::diagRowMul(C->__elem, A->__elem, std::min(*M, *K), *N);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 data_col = std::min(*K, *N);

      for(int r = 0; r < (int)*M; r++)
        uni10_elem_cast(&C->__elem[r*(*N)], &A->__elem[r*(*K)], data_col);

      uni10_linalg::diagColMul(C->__elem, B->__elem, *M, data_col);

    }
    else{

      uni10_uint64 min = std::min(A->__elemNum, B->__elemNum);
      uni10_elem_copy(C->__elem, A->__elem, min*sizeof(uni10_complex128));

      uni10_linalg::vectorMul(C->__elem, B->__elem, min);

    }

  }

  void matrixDot(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_linalg::matrixDot(A->__elem, B->__elem, *M, *N, *K, C->__elem);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_elem_cast(C->__elem, B->__elem, B->__elemNum);
      uni10_linalg::diagRowMul(C->__elem, A->__elem, std::min(*M, *K), *N);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 data_col = std::min(*K, *N);

      for(int r = 0; r < (int)*M; r++)
        uni10_elem_copy(&C->__elem[r*(*N)], &A->__elem[r*(*K)], data_col*sizeof(uni10_complex128));

      uni10_linalg::diagColMul(C->__elem, B->__elem, *M, data_col);

    }
    else{

      uni10_uint64 min = std::min(A->__elemNum, B->__elemNum);
      uni10_elem_copy(C->__elem, A->__elem, min*sizeof(uni10_complex128));

      uni10_linalg::vectorMul(C->__elem, B->__elem, min);

    }

  }

}
