#include "5pess_kagome.h"
#include "5pess_tools.h"

template<typename T>
PESS_5<T>::PESS_5(const UniTensor<T>& twoSiteH,  const pess5_paras& paras, const map<string, Network_dev*>& _net_list): UsNum(4), net_list(_net_list){

  fprintf(stdout, "\n");
  paras.print_info();
  fprintf(stdout, "\n");

  cout << twoSiteH;

  H_FOR_5PESS(twoSiteH);

  dim = H.bond(0).dim();
  D = paras.D;
  max_N   = paras.max_N;
  tau     = paras.tau;
  eps     = paras.eps;
  cut_off = paras.cut_off;
  measure_per_n_iter = paras.measure_per_n_iter;

  max_D = (uni10_int)cut_off == -1 ? D : paras.max_D;
  this->init();

}

template<typename T>
void PESS_5<T>::H_FOR_5PESS(const UniTensor<T>& twoSiteH){

  uni10_int physdim = twoSiteH.bond(0).dim();
  vector<Bond> bondsH(10, Bond(BD_IN, physdim));
  for(uni10_int i = 5; i < 10; i++)
    bondsH[i] = Bond(BD_OUT, physdim);

  H.assign(bondsH);
  vector<uni10_int> ori_labels = H.label();
  Matrix<T>  id(physdim, physdim, true);
  id.identity();
  Matrix<T>  HELEM = otimes( otimes( otimes(twoSiteH.getBlock(), id), id), id);
  H.putBlock(HELEM);

  int per1_labels[] = {2, 0, 1, 3, 4, 7, 5, 6, 8, 9};
  permute(H,per1_labels, 5,INPLACE);
  HELEM += H.getBlock();

  int per2_labels[] = {2, 3, 0, 1, 4, 7, 8, 5, 6, 9};
  permute(H,per2_labels, 5,INPLACE);
  HELEM += H.getBlock();

  int per3_labels[] = {2, 3, 4, 0, 1, 7, 8, 9, 5, 6};
  permute(H,per3_labels, 5,INPLACE);
  HELEM += H.getBlock();

  int per4_labels[] = {0, 2, 1, 3, 4, 5, 7, 6, 8, 9};
  permute(H,per4_labels, 5,INPLACE);
  HELEM += H.getBlock();

  int per5_labels[] = {2, 3, 0, 4, 1, 7, 8, 5, 9, 6};
  permute(H,per5_labels, 5,INPLACE);
  HELEM += H.getBlock();

  H.putBlock(HELEM);
  H.setLabel(ori_labels);

}

template<typename T>
void PESS_5<T>::init(){

  vector<Bond> U_bonds(3, Bond(BD_IN, D));
  U_bonds[0] = H.bond(0);

  Us.assign(UsNum, UniTensor<T> (U_bonds));

  for(uni10_uint64 i = 0; i < Us.size(); i++)
    Us[i].randomize();

  U_bonds.push_back(U_bonds[1]);
  U_bonds.push_back(U_bonds[1]);
  Cs.assign(2, UniTensor<T> (U_bonds));

  for(uni10_uint64 i = 0; i < Cs.size(); i++)
    Cs[i].randomize();

  Matrix<T> _unitLs(D, D, true);
  unitLs.assign(2, _unitLs);
  uni10_rand(_unitLs, uni10_mt19937, uni10_normal, 0, 1, 777);

  for(uni10_uint64 i = 0; i < unitLs.size(); i++)
    unitLs[i] = _unitLs;

  Ls.assign(UsNum, unitLs);

}

template<typename T>
PESS_5<T>::~PESS_5(){

}

template<typename T>
void PESS_5<T>::setHamiltonian(const UniTensor<T>& twoSiteH){

  H_FOR_5PESS(twoSiteH);

}

template<typename T>
UniTensor<T> PESS_5<T>::get_gate(){

  UniTensor<T>  expH(H.bond());
  expH.putBlock(exph(-tau, H.const_getBlock()));
  return expH;

}

template<typename T>
void PESS_5<T>::Optimize(){

  fprintf(stdout, "Updating the matrix product states: \n\n");
  progressbar(0, 0, max_N, true);

  uni10_double64 preGE = 0, GE , err; 

  UniTensor<T> expH = get_gate();

  for(uni10_uint64 i = 0; i < max_N; i++){

    update_driver(0, expH);
    update_driver(1, expH);

    err = (preGE - GE)/6;

    if(i%measure_per_n_iter==0){ 

      GE = 0;
      update_driver(0, expH);
      GE += measure_driver(0, H);
      update_driver(1, expH);
      GE += measure_driver(1, H);
      err = (preGE - GE)/6;

      progressbar(i+measure_per_n_iter, 0, max_N);
      cout.precision(8);
      cout.setf(ios::fixed, ios::floatfield);
      cout <<  ", ge: " << GE/6. << "  err: " << err/measure_per_n_iter <<"\r";
      std::cout.flush(); 
      preGE = GE;

    }  

  }

  fprintf(stdout, "\n\n");

}

template<typename T>
void PESS_5<T>::update_driver( int dir, const UniTensor<T>& expH){

  bondscat(dir, Us, Ls);
  permuteUs(dir, Us);
  vector<uni10_int> ori_labels = Us[0].label();

  std::map<string, Network_dev*>::iterator it = net_list.find("theta");

  UniTensor<T>  theta;
  contract_args(theta, *it->second, Us[0], Us[1], Us[2], Us[3], Cs[dir], expH);

  vector< uni10_int > group_labels = theta.label();
  vector< uni10_int > groups(4,2);
  vector< Matrix<T> > svdLs;
  vector< UniTensor<T> > svdUs;
  UniTensor< T > Core;
  hosvd(theta, group_labels, groups, svdUs, Core, svdLs, INPLACE);
  svdUs.push_back(Core);

  //truncate Core
  truncateLUs(dir, D, Us, Ls, svdLs, svdUs, cut_off);
  UniTensor<T>  trunC = theta;

  for(uni10_uint64 i = 0; i < svdUs.size()-1; i++){
    trunC = contract(trunC ,svdUs[i], true);
  }

  uni10_double64 nrm = norm(trunC.getBlock());

  UniTensor<T> pT;
  pT = trunC * (1.0 / nrm);
  permute(Cs[dir], pT, trunC.bondNum(), INPLACE);

  bondsrm(dir, Us, Ls);

}

template<typename T> template<typename U>
uni10_double64 PESS_5<T>::measure_driver(int dir, const UniTensor<U> & Ob){

  vector<UniTensor<T> > _Us = Us;

  bondscat(dir, _Us, Ls);
  permuteUs(dir, _Us);

  UniTensor<T>  S, SOST;

  std::map<string, Network_dev*>::iterator it0 = net_list.find("state");

  it0->second->putTensor("Core", Cs[dir]);
  for(uni10_int i = 0; i < _Us.size(); i++)
    it0->second->putTensor(i+1, _Us[i]); // C
  it0->second->launch(S);

  std::map<string, Network_dev*>::iterator it1 = net_list.find("measure");
  it1->second->putTensor("S", S);
  it1->second->putTensor("ST", S);
  it1->second->putTensor("Ob", Ob);
  it1->second->launch(SOST);

  return GETREAL(SOST[0]) / norm(S.getBlock());

}

template class PESS_5<uni10_double64>;
template class PESS_5<uni10_complex128>;
