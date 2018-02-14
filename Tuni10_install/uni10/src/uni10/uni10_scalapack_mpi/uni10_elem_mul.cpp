#include "uni10/uni10_scalapack_mpi/uni10_elem_linalg_scalapack_mpi.h"

namespace uni10{

  void matrixMul(uni10_elem_double64* A, uni10_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* N){

    if(isAdiag && !isBdiag){
      for(uni10_uint64 i = 0; i < A->__elemNum; i++)
        A->__elem[i] *= B->__elem[i*(*N)+i];
    }
    else if(!isAdiag && isBdiag){

      *isAdiag = true;

      A->__elemNum = B->__elemNum;

      uni10_double64* _elem = (uni10_double64*)malloc(B->__elemNum * sizeof(uni10_double64));

      for(uni10_uint64 i = 0; i < B->__elemNum; i++)
        _elem[i] = A->__elem[i*(*N)+i] * B->__elem[i];

      if(A->__elem != NULL)
        uni10_elem_free(A->__elem, A->__elemNum*sizeof(uni10_double64));

      A->__elem = _elem;
    }
    else
      uni10_linalg::vectorMul(A->__elem, B->__elem, B->__elemNum);

  }

  void matrixMul(uni10_elem_complex128* A, uni10_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* N){

    if(isAdiag && !isBdiag){
      for(uni10_uint64 i = 0; i < A->__elemNum; i++)
        A->__elem[i] *= B->__elem[i*(*N)+i];
    }
    else if(!isAdiag && isBdiag){

      *isAdiag = true;

      A->__elemNum = B->__elemNum;

      uni10_complex128* _elem = (uni10_complex128*)malloc(B->__elemNum * sizeof(uni10_complex128));

      for(uni10_uint64 i = 0; i < B->__elemNum; i++)
        _elem[i] = A->__elem[i*(*N)+i] * B->__elem[i];

      if(A->__elem != NULL)
        uni10_elem_free(A->__elem, A->__elemNum*sizeof(uni10_complex128));

      A->__elem = _elem;
    }
    else
      uni10_linalg::vectorMul(A->__elem, B->__elem, B->__elemNum);


  }

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* isAdiag, const uni10_elem_double64* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C){

    if( *isAdiag && !*isBdiag ){

      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++){
        C->__elem[i] = A->__elem[i] * B->__elem[i * (*N) + i];
      }

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i] = B->__elem[i] * A->__elem[i * (*N) + i];

    }
    else{

      uni10_elem_copy(C->__elem, A->__elem,  A->__elemNum*sizeof(uni10_double64));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }

  }

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* isAdiag, const uni10_elem_complex128* B, uni10_const_bool* isBdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C){

    if( *isAdiag && !*isBdiag ){

      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++){
        C->__elem[i] = A->__elem[i] * B->__elem[i * (*N) + i];
      }

    }
    else if( !*isAdiag && *isBdiag ){

      uni10_uint64 min = std::min(*M, *N);
      for(int i = 0; i < (int)min; i++)
        C->__elem[i] = B->__elem[i] * A->__elem[i * (*N) + i];

    }
    else{

      uni10_elem_copy(C->__elem, A->__elem,  A->__elemNum*sizeof(uni10_complex128));
      uni10_linalg::vectorMul(C->__elem, B->__elem, C->__elemNum);

    }

  }

}
