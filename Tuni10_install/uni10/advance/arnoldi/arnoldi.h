#ifndef __UNI10_ARNOLDI_H__
#define __UNI10_ARNOLDI_H__

#include "uni10/uni10.hpp"

#include "customized_structure.h"

using namespace std;
using namespace uni10;

template<typename T, typename U>
class uni10_arnoldi_paras{

  public:

    uni10_arnoldi_paras(): paras(NULL), customID(NULL), Mout(NULL){}

    uni10_arnoldi_paras(U& _paras){

      this->paras = &_paras;
      this->customID = &paras->customID;
      this->MW = &paras->MW;
      this->Mout = paras->Mout;

    }

    ~uni10_arnoldi_paras(){

    }

    void setParas(U& _para){

      this->paras = &_para;
      this->customID = &paras->customID;
      this->MW = &paras->MW;
      this->Mout = paras->Mout;

    };

    void arnoldi_mul(){
      paras->arnoldi_mul();
    }

    template<typename _T, typename _U>
      friend void Arnoldi(std::vector<uni10_complex128>& eigvs, const uni10_uint64 eigvNum, uni10_arnoldi_paras<_T, _U>& _paras, uni10_uint64& n, 
          uni10_uint64 maxN, uni10_uint64 minN, 
          uni10_int32& info, const uni10_double64 cut_off);

  private:
    U* paras; 
    uni10_int32* customID;
    uni10::Matrix<T>* MW;
    uni10::Matrix<T>* Mout;

};

template<typename T, typename U>
uni10_complex128 Arnoldi(uni10_arnoldi_paras<T, U>& _paras, uni10_uint64& n, uni10_uint64 maxN, uni10_uint64 minN,
    uni10_int32& info, const uni10_double64 cut_off=1E-14);

template<typename T, typename U>
void Arnoldi(std::vector<uni10_complex128>& eigvs, const uni10_uint64 eigvNum, uni10_arnoldi_paras<T, U>& _paras, uni10_uint64& n, 
    uni10_uint64 maxN, uni10_uint64 minN, 
    uni10_int32& info, const uni10_double64 cut_off=1E-14);

template<typename T, typename U>
void Arnoldi(uni10::Matrix<uni10_complex128>& eigvs, const uni10_uint64 eigvNum, uni10_arnoldi_paras<T, U>& _paras, uni10_uint64& n, 
    uni10_uint64 maxN, uni10_uint64 minN, 
    uni10_int32& info, const uni10_double64 cut_off=1E-14);

template<typename T, typename U>
uni10_complex128 Arnoldi(uni10_arnoldi_paras<T, U>& _paras, uni10_uint64& n,  uni10_uint64 maxN, uni10_uint64 minN, uni10_int32& info, const uni10_double64 cut_off){

  std::vector<uni10_complex128> eigvs;
  Arnoldi(eigvs, 1, _paras, n, maxN, minN, info, cut_off);
  return eigvs[0];

}

template<typename T, typename U>
void Arnoldi(uni10::Matrix<uni10_complex128>& eigvs, const uni10_uint64 eigvNum, uni10_arnoldi_paras<T, U>& _paras, uni10_uint64& n, 
    uni10_uint64 maxN, uni10_uint64 minN, 
    uni10_int32& info, const uni10_double64 cut_off){

  std::vector<uni10_complex128> _eigvs;
  Arnoldi(_eigvs, eigvNum, _paras, n, maxN, minN, info, cut_off);
  eigvs.assign(eigvNum, eigvNum, true);
  eigvs.setElem(_eigvs);

}

template<typename T, typename U>
void Arnoldi(std::vector<uni10_complex128>& eigvs, const uni10_uint64 eigvNum, uni10_arnoldi_paras<T, U>& _paras, uni10_uint64& n, 
    uni10_uint64 maxN, uni10_uint64 minN, 
    uni10_int32& info, const uni10_double64 cut_off){

  info = 0;  // Checking the convergence.
  uni10::Matrix<T> V  = *_paras.Mout;
  _paras.arnoldi_mul();
  do{
    n++;
    for(uni10_uint64 i = 0; i < n; i++){
        
    } 
    dot(V, , V)
  }while(n < maxN && info != -1)


}


#endif
