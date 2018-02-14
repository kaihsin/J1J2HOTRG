#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"

namespace uni10{
  // A += B
  void matrixAdd(uni10_elem_double64* A, uni10_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N){

    if(isAdiag && !isBdiag){

      uni10_error_msg(true, "%s", "Developing");

    }
    else if(!isAdiag && isBdiag){

      uni10_error_msg(true, "%s", "Developing");

    }
    else
      uni10_linalg::vectorAdd(A->__elem, B->__elem, B->blockrow*B->blockcol);

  }

  void matrixAdd(uni10_elem_complex128* A, uni10_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N){

    if(isAdiag && !isBdiag){

      uni10_error_msg(true, "%s", "Developing");

    }
    else if(!isAdiag && isBdiag){

      uni10_error_msg(true, "%s", "Developing");

    }
    else
      uni10_linalg::vectorAdd(A->__elem, B->__elem, B->blockrow*B->blockcol);

  }

  void matrixAdd(uni10_elem_complex128* A, uni10_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N){

    if(isAdiag && !isBdiag){

      uni10_error_msg(true, "%s", "Developing");

    }
    else if(!isAdiag && isBdiag){

      uni10_error_msg(true, "%s", "Developing");

    }
    else{

      uni10_linalg::vectorAdd(A->__elem, B->__elem, B->blockrow*B->blockcol);

    }

  }

  // C = A + B
  void matrixAdd(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy(C->__elem, B->__elem, B->blockrow*B->blockcol * sizeof(uni10_double64) );
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->blockrow*C->blockcol);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_error_msg(true, "%s", "Developing");

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_error_msg(true, "%s", "Developing");

    }
    else{

      uni10_elem_copy(C->__elem, B->__elem, B->blockrow*B->blockcol * sizeof(uni10_double64));
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->blockrow*C->blockcol);

    }

  }

  void matrixAdd(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128));
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->blockrow*C->blockcol);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_error_msg(true, "%s", "Developing");

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_error_msg(true, "%s", "Developing");

    }
    else{

      uni10_elem_copy(C->__elem, B->__elem,  C->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->blockrow*C->blockcol);

    }

  }

  void matrixAdd(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy(C->__elem, B->__elem, B->__elemNum * sizeof(uni10_complex128) );
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->blockrow*C->blockcol);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_error_msg(true, "%s", "Developing");

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_error_msg(true, "%s", "Developing");

    }
    else{

      uni10_elem_copy(C->__elem, B->__elem,  C->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->blockrow*C->blockcol);

    }

  }

  void matrixAdd(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( !*isAdiag && !*isBdiag ){

      uni10_elem_copy(C->__elem, A->__elem, A->__elemNum * sizeof(uni10_complex128) );
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->blockrow*C->blockcol);

    }
    else if( *isAdiag && !*isBdiag ){

      uni10_error_msg(true, "%s", "Developing");

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_error_msg(true, "%s", "Developing");

    }
    else{

      uni10_elem_copy(C->__elem, A->__elem, A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorAdd(C->__elem, A->__elem, C->blockrow*C->blockcol);

    }

  }

}
