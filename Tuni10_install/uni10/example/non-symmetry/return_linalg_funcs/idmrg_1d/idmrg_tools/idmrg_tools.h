#ifndef __IDMRG_TOOLS_H__
#define __IDMRG_TOOLS_H__

#include "../../../../../advance/lanczos/lanczos.h"

// Customized lanczos mutiplication.
//
uni10::UniTensor<double> mpoMatVec(
    std::vector<uni10::UniTensor<double>>& mpoH,
    uni10::UniTensor<double>& psi );

template<class T>
struct uni10_lanczos_custom3_paras{

  uni10_lanczos_custom3_paras(
      std::vector<uni10::UniTensor<T> >& _mpoH, uni10::UniTensor<T>& psi,
      uni10::Matrix<T>& psi_blk ): mpoH(_mpoH) {
    customID = 3;
    Mout = &psi_blk;
    psi_bd = psi.bond();
  }

  // Destructor
  ~uni10_lanczos_custom3_paras(){};

  // Customized multiplication.
  void lanczos_mul(){
    uni10::UniTensor<T> psi(psi_bd);
    psi.putBlock(*Mout);
    uni10::UniTensor<T> MV = mpoMatVec(mpoH, psi);
    MW = MV.getBlock();
  };

  std::vector<uni10::UniTensor<T> > mpoH;
  std::vector<uni10::Bond> psi_bd;

  uni10_int32 customID;
  uni10::Matrix<T> MW;
  uni10::Matrix<T>* Mout;
};

template struct uni10_lanczos_custom3_paras<double>;



uni10_uint64 LanczosEigh(
    std::vector<uni10::UniTensor<double>>& mpoH,
    uni10::UniTensor<double>& psi,
    double& E0, int max_iter, double err_tol );

uni10_uint64 LanczosEigh(
    uni10::UniTensor<double>& op,
    uni10::UniTensor<double>& psi,
    double& E0, int max_iter, double err_tol );

uni10::UniTensor<double> netLG( uni10::UniTensor<double> la, uni10::UniTensor<double> ga );
uni10::UniTensor<double> netGL( uni10::UniTensor<double> ga, uni10::UniTensor<double> la );
uni10::UniTensor<double> netLGLGL(
    uni10::UniTensor<double> ll, uni10::UniTensor<double> ga,
    uni10::UniTensor<double> la, uni10::UniTensor<double> gb, uni10::UniTensor<double> lb );

uni10::UniTensor<double> theta( uni10::UniTensor<double> ket, uni10::UniTensor<double> op );
uni10::UniTensor<double> expVal( uni10::UniTensor<double> ket, uni10::UniTensor<double> op );
uni10::UniTensor<double> tenInv( uni10::UniTensor<double>& ten );

uni10::UniTensor<double> contrMPOLR( std::vector<uni10::UniTensor<double>>& mpo );

void renormMPOL(
    std::vector<uni10::UniTensor<double>>& mpoH,
    uni10::UniTensor<double> ket, bool initial = false );

void renormMPOR(
    std::vector<uni10::UniTensor<double>>& mpoH,
    uni10::UniTensor<double> ket, bool initial = false );


#endif
