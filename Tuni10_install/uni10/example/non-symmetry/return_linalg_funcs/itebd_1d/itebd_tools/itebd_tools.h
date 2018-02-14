#ifndef __ITEBD_TOOLS_H__
#define __ITEBD_TOOLS_H__

#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

template<typename T>
void bondcat(UniTensor<T>& Tout, const Matrix<T>& L, uni10_int32 bidx);

template<typename T>
void bondrm(UniTensor<T>& Tout, const Matrix<T>& L, uni10_int32 bidx);

template<typename T>
void bondcat(UniTensor<T>& Tout, const Matrix<T>& L, uni10_int32 bidx){

  uni10_int32 inBondNum = Tout.inBondNum();
  vector<uni10_int32> labels = Tout.label();
  vector<uni10_int32> per_labels = labels;
  uni10_int32 l = labels[bidx];
  per_labels.erase(per_labels.begin() + bidx);
  per_labels.insert(per_labels.begin(), l);

  UniTensor<T> T_c = permute(Tout, per_labels, 1);
  T_c.putBlock((dot(L, T_c.getBlock())));
  Tout = permute( T_c, labels, inBondNum);
}

template<typename T>
void bondrm(UniTensor<T>& Tout, const Matrix<T>& L, uni10_int32 bidx){

  Matrix<T> invL = L;
  for(uni10_uint64 i=0; i!=L.col(); i++){
    invL[i] = invL[i] == 0.0 ? 0.0 : ( 1.0 / invL[i]);
  }
  bondcat(Tout, invL, bidx);

}

#endif
