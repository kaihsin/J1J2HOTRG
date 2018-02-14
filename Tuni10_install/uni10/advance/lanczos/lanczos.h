#ifndef __UNI10_LANCZOS_H__
#define __UNI10_LANCZOS_H__

#include "uni10/uni10.hpp"

#include "customized_structure.h"

template<typename T, typename U>
class uni10_lanczos_paras{

  public:

    uni10_lanczos_paras(): paras(NULL), customID(NULL), Mout(NULL){}

    uni10_lanczos_paras(U& _paras){

      this->paras = &_paras;
      this->customID = &paras->customID;
      this->MW = &paras->MW;
      this->Mout = paras->Mout;

    }

    ~uni10_lanczos_paras(){

    }

    void setParas(U& _para){

      this->paras = &_para;
      this->customID = &paras->customID;
      this->MW = &paras->MW;
      this->Mout = paras->Mout;

    };

    void lanczos_mul(){
      paras->lanczos_mul();
    }

    /*
    template<typename _T, typename _U>
      friend _T Lanczos(uni10_lanczos_paras<_T, _U>& paras, uni10_uint64& n, uni10_uint64 maxN, 
          uni10_uint64 minN, uni10_int32& info, const uni10_double64 cut_off);
    */
    template<typename _T, typename _U>
      friend void Lanczos(std::vector<uni10_double64>& eigvs, const uni10_uint64 eigvNum, uni10_lanczos_paras<_T, _U>& _paras, uni10_uint64& n, 
          uni10_uint64 maxN, uni10_uint64 minN, 
          uni10_int32& info, const uni10_double64 cut_off);

  private:
    U* paras; 
    uni10_int32* customID;
    uni10::Matrix<T>* MW;
    uni10::Matrix<T>* Mout;

};

template<typename T, typename U>
T Lanczos(uni10_lanczos_paras<T, U>& _paras, uni10_uint64& n, uni10_uint64 maxN, uni10_uint64 minN, 
    uni10_int32& info, const uni10_double64 cut_off=1E-12);


template<typename T, typename U>
void Lanczos(std::vector<uni10_double64>& eigvs, const uni10_uint64 eigvNum, uni10_lanczos_paras<T, U>& _paras, uni10_uint64& n, 
    uni10_uint64 maxN, uni10_uint64 minN, 
    uni10_int32& info, const uni10_double64 cut_off=1E-12);

template<typename T, typename U>
void Lanczos(uni10::Matrix<uni10_double64>& eigvs, const uni10_uint64 eigvNum, uni10_lanczos_paras<T, U>& _paras, uni10_uint64& n, 
    uni10_uint64 maxN, uni10_uint64 minN, 
    uni10_int32& info, const uni10_double64 cut_off=1E-12);

template<typename T, typename U>
T Lanczos(uni10_lanczos_paras<T, U>& _paras, uni10_uint64& n,  uni10_uint64 maxN, uni10_uint64 minN, uni10_int32& info, const uni10_double64 cut_off){

  std::vector<uni10_double64> eigvs;
  Lanczos(eigvs, 1, _paras, n, maxN, minN, info, cut_off);
  return eigvs[0];

}

template<typename T, typename U>
void Lanczos(uni10::Matrix<uni10_double64>& eigvs, const uni10_uint64 eigvNum, uni10_lanczos_paras<T, U>& _paras, uni10_uint64& n, 
    uni10_uint64 maxN, uni10_uint64 minN, 
    uni10_int32& info, const uni10_double64 cut_off){

  std::vector<uni10_double64> _eigvs;
  Lanczos(_eigvs, eigvNum, _paras, n, maxN, minN, info, cut_off);
  eigvs.assign(eigvNum, eigvNum, true);
  eigvs.setElem(_eigvs);

}

template<typename T, typename U>
void Lanczos(std::vector<uni10_double64>& eigvs, const uni10_uint64 eigvNum, uni10_lanczos_paras<T, U>& _paras, uni10_uint64& n, 
    uni10_uint64 maxN, uni10_uint64 minN, 
    uni10_int32& info, const uni10_double64 cut_off){


  info = 0;  // Checking the convergence.

  std::vector<T> Es;

  uni10::Matrix<T>& V = *_paras.Mout;

  uni10_uint64 Rnum = V.row();
  uni10_uint64 Cnum = V.col();

  // uni10_error_msg(maxN>Rnum, "%s", "The maximum iteratoin number can't be larger than the row number of V");
  // uni10_error_msg(minN>Rnum, "%s", "The minimum iteratoin number can't be larger than the row number of V");

  uni10::Matrix<T> pre_v(Rnum, Cnum); 

  uni10::Matrix<uni10_double64> vs(1, Rnum); 
  vs.elem_enforce().copy(0, V.const_elem_enforce(), Rnum);

  T alpha = 0.0; T beta = 0.0;

  uni10::UELEM(uni10_elem, _package, _type)<T> Alphas(1, maxN+1); 
  uni10::UELEM(uni10_elem, _package, _type)<T> Betas(1, maxN); 
  uni10::UELEM(uni10_elem, _package, _type)<T> As(1, maxN+1);
  uni10::UELEM(uni10_elem, _package, _type)<T> Bs(1, maxN); 

  uni10_double64 err = 1.0;

  do{

    _paras.lanczos_mul();
    uni10::Matrix<T> W = *_paras.MW; 
    uni10::Matrix<T> WT;
    transpose(WT, W, uni10::INPLACE); 
    uni10::Matrix<T> WTV;
    dot(WTV, WT, V, uni10::INPLACE);
    
    alpha = WTV[0];
    Alphas[n] = alpha;
    // Add function add_args.
    W = W + (-1)*alpha*V + (-1)*beta*pre_v;
    beta = norm(W); 
    Betas[n] = beta;

    n++; info=n;
    pre_v = V;
    V = (1.0/beta) * W; 

    As.copy(0, Alphas, n);
    Bs.copy(0, Betas , n-1);

    trimatrixEigh(&As, &Bs, &n);
    Es.push_back(As[0]);

    if(n >= minN)
      err = fabs((As[0]-Es[n-minN])/(fabs(As[0]) > 1. ? fabs(As[0]) : 1));
    
    if(err < cut_off && n >= minN){
      info = -1;
      break;
    }

    if(n < maxN){
      vs.row_enforce() += 1;
      vs.elem_enforce().catElem(V.const_elem_enforce());
    }

  }while(n < maxN);

  As.copy(0, Alphas, n);
  Bs.copy(0, Betas, n-1);

  uni10::Matrix<uni10_double64> mus(n, n); 
  trimatrixEigh(&As, &Bs, &n, &mus.elem_enforce(), &n);
  resize(mus, eigvNum, n, uni10::INPLACE);
  dot(V, mus, vs, uni10::INPLACE); 
  uni10::Matrix<uni10_double64> vs_sub(1, Rnum);
 
  for(uni10_uint64 i = 0; i < eigvNum; i++){
    vs_sub.elem_enforce().copy(0, V.const_elem_enforce(), i*Rnum, Rnum);
    vs_sub *= 1.0/ norm(vs_sub);
    V.elem_enforce().copy(i*Rnum, vs_sub.const_elem_enforce(), Rnum);
  }

  transpose(V, uni10::INPLACE);
  eigvs.assign(&As[0], &As[0+eigvNum]); 

  if(info != -1)
    fprintf(stdout, " not converge \n");

}


#endif
