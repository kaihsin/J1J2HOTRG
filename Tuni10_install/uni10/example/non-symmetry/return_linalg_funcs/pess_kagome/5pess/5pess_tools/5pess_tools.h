#ifndef __ITEBD_TOOLS_H__
#define __ITEBD_TOOLS_H__

#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

uni10_double64 GETREAL(uni10_double64 R);

uni10_double64 GETREAL(uni10_complex128 C);

template<typename T>
void bondcat(UniTensor<T>& Tout, const Matrix<T>& L, uni10_int32 bidx);

template<typename T>
void bondrm(UniTensor<T>& Tout, const Matrix<T>& L, uni10_int32 bidx);

template<typename T>
void bondscat(uni10_int dir, vector<UniTensor<T> >& Us, const vector< vector<Matrix<T> > >& Ls);

template<typename T>
void bondsrm(uni10_int dir, vector<UniTensor<T> >& Us, const vector< vector<Matrix<T> > >& Ls);

template<typename T>
void permuteUs(uni10_int dir, vector<UniTensor<T> >& Us);

template<typename T>
int truncateLUs(uni10_int dir, uni10_int chi, vector<UniTensor<T> >& Us, vector< vector<Matrix<T> > >& Ls, vector<Matrix<T> >& svdLs, 
    vector<UniTensor<T> >& svdUs, const uni10_double64& accuracy);



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


template<typename T>
void bondscat(uni10_int dir, vector<UniTensor<T> >& Us, const vector< vector<Matrix<T> > >& Ls){

  for(uni10_uint64 i = 0; i < Us.size(); i++)
    for(uni10_uint64 j = 0; j < Ls[i].size(); j++)
      if( j != dir )
        bondcat(Us[i], Ls[i][j], j+1);

}

template<typename T>
void bondsrm(uni10_int dir, vector<UniTensor<T> >& Us, const vector< vector<Matrix<T> > >& Ls){

  for(uni10_uint64 i = 0; i < Us.size(); i++)
    for(uni10_uint64 j = 0; j < Ls[i].size(); j++)
      if( j != dir )
        bondrm(Us[i], Ls[i][j], j+1);

}

template<typename T>
void permuteUs(uni10_int dir, vector<UniTensor<T> >& Us){

  if(dir == 0){
    uni10_int per_labels[] = {0, 2, 1};
    for(uni10_int i = 0; i < Us.size(); i++)
      permute(Us[i], per_labels, 2, INPLACE);
  }
  if(dir == 1){
    uni10_int per_labels[] = {0, 1, 2};
    for(uni10_int i = 0; i < Us.size(); i++)
      permute(Us[i], per_labels, 2, INPLACE);
  }

}

template<typename T>
int truncateLUs(uni10_int dir, uni10_int D, vector<UniTensor<T> >& Us, vector< vector<Matrix<T> > >& Ls, vector<Matrix<T> >& svdLs, 
    vector<UniTensor<T> >& svdUs, const uni10_double64& cut_off){

  uni10_int deafult_labels[] = {0, 1, 2};

  vector<uni10_int> ori_labels = Us[0].label();
  vector<Bond> new_bonds;

  uni10_uint64 max_D = 0;

  for(uni10_uint64 i = 0; i < svdLs.size(); i++){

    new_bonds = Us[i].bond();

    uni10_double64 err = 1;
    uni10_int D_cut    = 0;

    vector<uni10_int> ori_svdUs_labels = svdUs[i].label();

    uni10_double64 Norm = norm(svdLs[i]);
    svdLs[i] *= (1.0/Norm);

    if(cut_off != -1){

      for(uni10_int j = 0; j < svdLs[i].col(); j++){

        err -= pow(GETREAL(svdLs[i][j]), 2);

        if(err < cut_off){
          D_cut = j + 1;
          break;
        }

      }

      D_cut = D_cut <= D ? D : D_cut;

    }else{

      D_cut = D;

    }

    max_D = max_D < D ? D : max_D;

    resize(Ls[i][dir],svdLs[i], D_cut, D_cut,INPLACE);
    new_bonds[2] = Bond(BD_OUT, D_cut);

    Matrix< T > blk;
    resize(blk, svdUs[i].getBlock(), svdUs[i].getBlock().row(), D_cut, INPLACE);
    svdUs[i].assign(new_bonds);
    svdUs[i].putBlock(blk);
    Us[i] = svdUs[i];
    svdUs[i].setLabel(ori_svdUs_labels);
    Us[i].setLabel(ori_labels);

    permute(Us[i], deafult_labels, 3, INPLACE);

  }

  return max_D;

}

#endif
