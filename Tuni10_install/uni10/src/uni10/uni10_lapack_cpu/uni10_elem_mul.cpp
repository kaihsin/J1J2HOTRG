#include "uni10/uni10_lapack_cpu/uni10_elem_linalg_lapack_cpu.h"

namespace uni10{

  void matrixMul(uni10_elem_double64* A, uni10_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64*M, const uni10_uint64* N){

    if(*isAdiag && !*isBdiag){

      uni10_linalg::matrix_diag_dense_mul(A->__elem, B->__elem, *M, *N);

    }
    else if(!*isAdiag && *isBdiag){

      uni10_elem_double64 MA(*A);
      A->init(*M, *N, true);
      uni10_linalg::matrix_dense_diag_mul(MA.__elem, B->__elem, *M, *N, A->__elem);

    }
    else
      uni10_linalg::vectorMul(A->__elem, B->__elem, B->__elemNum);

  }

  void matrixMul(uni10_elem_complex128* A, uni10_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64 *M, const uni10_uint64* N){

    if(*isAdiag && !*isBdiag){

      uni10_linalg::matrix_diag_dense_mul(A->__elem, B->__elem, *M, *N);

    }
    else if(!*isAdiag && *isBdiag){

      uni10_elem_complex128 MA(*A);
      A->init(*M, *N, true);
      uni10_linalg::matrix_dense_diag_mul(MA.__elem, B->__elem, *M, *N, A->__elem);

    }
    else
      uni10_linalg::vectorMul(A->__elem, B->__elem, B->__elemNum);


  }

  void matrixMul(uni10_elem_complex128* A, uni10_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64 *M, const uni10_uint64* N){

    if(*isAdiag && !*isBdiag){

      uni10_linalg::matrix_diag_dense_mul(A->__elem, B->__elem, *M, *N);

    }
    else if(!*isAdiag && *isBdiag){

      uni10_elem_complex128 MA(*A);
      A->init(*M, *N, true);
      uni10_linalg::matrix_dense_diag_mul(MA.__elem, B->__elem, *M, *N, A->__elem);

    }
    else
      uni10_linalg::vectorMul(A->__elem, B->__elem, B->__elemNum);


  }

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C){

    if( *isAdiag && !*isBdiag ){

      uni10_linalg::matrix_diag_dense_mul(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_linalg::matrix_dense_diag_mul(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else{

      uni10_elem_copy(C->__elem, A->__elem,  A->__elemNum*sizeof(uni10_double64));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( *isAdiag && !*isBdiag ){

      uni10_linalg::matrix_diag_dense_mul(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_linalg::matrix_dense_diag_mul(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else{

      uni10_elem_copy(C->__elem, A->__elem,  A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( *isAdiag && !*isBdiag ){

      uni10_linalg::matrix_diag_dense_mul(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_linalg::matrix_dense_diag_mul(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else{

      uni10_elem_copy(C->__elem, B->__elem,  B->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorMul(C->__elem, A->__elem, C->__elemNum);

    }

  }

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( *isAdiag && !*isBdiag ){

      uni10_linalg::matrix_diag_dense_mul(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_linalg::matrix_dense_diag_mul(A->__elem, B->__elem, *M, *N, C->__elem);

    }
    else{

      uni10_elem_copy(C->__elem, A->__elem,  A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }

  }

}
