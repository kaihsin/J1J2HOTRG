#ifndef __UNI10_HIGH_RANK_LINALG_OTIMES_H__
#define __UNI10_HIGH_RANK_LINALG_OTIMES_H__

#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_otimes.h"
namespace uni10{

  template<typename T>
    UniTensor<uni10_complex128> otimes(const UniTensor<uni10_complex128>& Ta, const UniTensor<T>& Tb);

  template<typename T>
    UniTensor<T> otimes(const UniTensor<uni10_double64>& Ta, const UniTensor<T>& Tb);

  template<typename T>
    Matrix<uni10_complex128> otimes(const Block<uni10_complex128>& Ta, const Block<T>& Tb);

  template<typename T>
    Matrix<T> otimes(const Block<uni10_double64>& Ta, const Block<T>& Tb);


  template<typename T>
    UniTensor<uni10_complex128> otimes(const UniTensor<uni10_complex128>& Ta, const UniTensor<T>& Tb){

      UniTensor<uni10_complex128> Tout;
      otimes(Tout, Ta, Tb, INPLACE);
      return Tout;

    }

  template<typename T>
    UniTensor<T> otimes(const UniTensor<uni10_double64>& Ta, const UniTensor<T>& Tb){

      UniTensor<T> Tout;
      otimes(Tout, Ta, Tb, INPLACE);
      return Tout;

    }

  template<typename T>
    Matrix<uni10_complex128> otimes(const Block<uni10_complex128>& Ma, const Block<T>& Mb){

      Matrix<uni10_complex128> Mout;
      otimes(Mout, Ma, Mb, INPLACE);
      return Mout;

    }

  template<typename T>
    Matrix<T> otimes(const Block<uni10_double64>& Ma, const Block<T>& Mb){

      Matrix<T> Mout;
      otimes(Mout, Ma, Mb, INPLACE);
      return Mout;

    }

};

#endif
