#ifndef __UNI10_ARNOLDI_CUSTMUL_H__
#define __UNI10_ARNOLDI_CUSTMUL_H__

#include "uni10/uni10.hpp"

template<typename T>
struct uni10_arnoldi_default_paras{

  // Default constructor
  uni10_arnoldi_default_paras(){
    customID = 0;
  };

  uni10_arnoldi_default_paras(const uni10::Matrix<T>& _Min, uni10::Matrix<T>& _Mout){
    customID = 0;
    Min  = _Min;
    Mout = &_Mout;
  };

  ~uni10_arnoldi_default_paras(){};

  void arnoldi_mul(){

    //uni10::Matrix<T> Mtmp(Mout);
    dot(MW, Min, *Mout, uni10::INPLACE);

  }

  uni10_int32 customID;
  uni10::Matrix<T> Min;
  uni10::Matrix<T> MW;
  uni10::Matrix<T>* Mout;

};

template<typename T>
struct uni10_arnoldi_custom1_paras{

  // Default constructor
  uni10_arnoldi_custom1_paras(){
    customID = 1;
  };

  // Add your customized constructor in the below function prototype.
  uni10_arnoldi_custom1_paras(uni10::UniTensor<T>& _EnvL, uni10::UniTensor<T>& _WL, uni10::UniTensor<T>& _WR, uni10::UniTensor<T>& _EnvR, 
      uni10::Network<T>& _arnoldi_net, uni10::Matrix<T>& GS): EnvL(_EnvL), WL(_WL), WR(_WR), EnvR(_EnvR){
    customID = 1;
    Mout = &GS;
    Lanczos_net = &_arnoldi_net;
    std::vector<uni10::Bond> bondsT;
    bondsT.push_back(uni10::Bond(uni10::BD_IN, 1));
    bondsT.push_back(uni10::Bond(uni10::BD_IN, 2));
    bondsT.push_back(uni10::Bond(uni10::BD_IN, 2));
    bondsT.push_back(uni10::Bond(uni10::BD_IN, 1));
    Tout.assign(bondsT);

  };

  // Destructor
  ~uni10_arnoldi_custom1_paras(){};

  // Customized multiplication.
  void arnoldi_mul(){

    Tout.putBlock(*Mout);
    Lanczos_net->putTensor("EnvL", EnvL);
    Lanczos_net->putTensor("WL", WL);
    Lanczos_net->putTensor("WR", WR);
    Lanczos_net->putTensor("EnvR", EnvR);
    Lanczos_net->putTensor("GS", Tout);
    TW = Lanczos_net->launch();
    MW = TW.getBlock();


  }
  //  
  //
  // Add the parameters you need
  uni10::UniTensor<T> EnvL, WL, WR, EnvR, TW, Tout;
  uni10::Network<T>* Lanczos_net; 
  //
  //
  uni10_int32 customID;
  uni10::Matrix<T> MW;
  uni10::Matrix<T>* Mout;

};

template<typename T>
struct uni10_arnoldi_custom2_paras{

  uni10_arnoldi_custom2_paras(){
    customID = 2;

  }

  // Destructor
  ~uni10_arnoldi_custom2_paras(){};

  // Customized multiplication.
  void arnoldi_mul(){


  };
  //  
  //
  // Add the parameters you need
  // 
  //
  uni10_int32 customID;
  uni10::Matrix<T> MW;
  uni10::Matrix<T>* Mout;

};

#endif
