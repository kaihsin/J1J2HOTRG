#include "uni10/uni10_cusolver_gpu/uni10_elem_linalg_cusolver_gpu.h"

namespace uni10{

  void matrixSub(uni10_elem_double64* A, uni10_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N){

    if(*isAdiag && !*isBdiag){

      uni10_elem_double64 DA(*A);
      A->init(*M, *N, false);
      uni10_linalg::matrix_diag_dense_sub(DA.__elem, B->__elem, *M, *N, A->__elem);

    }
    else if(!*isAdiag && *isBdiag){

      uni10_linalg::matrix_dense_diag_sub(A->__elem, B->__elem, *M, *N);

    }
    else
      uni10_linalg::vectorSub(A->__elem, B->__elem, B->__elemNum);

  }

  void matrixSub(uni10_elem_complex128* A, uni10_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N){

    if(*isAdiag && !*isBdiag){

      uni10_elem_complex128 DA(*A);
      A->init(*M, *N, false);
      uni10_linalg::matrix_diag_dense_sub(DA.__elem, B->__elem, *M, *N, A->__elem);

    }
    else if(!*isAdiag && *isBdiag){

      uni10_linalg::matrix_dense_diag_sub(A->__elem, B->__elem, *M, *N);

    }
    else
      uni10_linalg::vectorSub(A->__elem, B->__elem, B->__elemNum);

  }

  void matrixSub(uni10_elem_complex128* A, uni10_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N){

    if(*isAdiag && !*isBdiag){

      uni10_elem_complex128 DA(*A);
      A->init(*M, *N, false);
      uni10_linalg::matrix_diag_dense_sub(DA.__elem, B->__elem, *M, *N, A->__elem);

    }
    else if(!*isAdiag && *isBdiag){

      uni10_linalg::matrix_dense_diag_sub(A->__elem, B->__elem, *M, *N);

    }
    else
      uni10_linalg::vectorSub(A->__elem, B->__elem, B->__elemNum);

  }

  void matrixSub(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C){

    if( *isAdiag && !*isBdiag ){

      uni10_linalg::matrix_diag_dense_sub(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_linalg::matrix_dense_diag_sub(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else{

      uni10_elem_copy(C->__elem, A->__elem,  C->__elemNum*sizeof(uni10_double64));
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixSub(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( *isAdiag && !*isBdiag ){

      uni10_linalg::matrix_diag_dense_sub(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_linalg::matrix_dense_diag_sub(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else{

      uni10_elem_copy(C->__elem, A->__elem,  A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixSub(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( *isAdiag && !*isBdiag ){

      uni10_linalg::matrix_diag_dense_sub(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_linalg::matrix_dense_diag_sub(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else{

      uni10_elem_cast(C->__elem, A->__elem, C->__elemNum);
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixSub(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( *isAdiag && !*isBdiag ){

      uni10_linalg::matrix_diag_dense_sub(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_linalg::matrix_dense_diag_sub(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else{

      uni10_elem_copy(C->__elem, A->__elem, A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorSub(C->__elem, B->__elem, C->__elemNum);

    }

  }

}
