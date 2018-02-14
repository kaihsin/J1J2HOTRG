#include "uni10/uni10.hpp"

#include "iostream"
#include "string"

using namespace std;
using namespace uni10;


int main(){

  int chi = 3;
  int D = 2;

  Bond bdi_chi(BD_IN, chi);
  Bond bdi_D(BD_IN, D);

  Bond bdo_chi(BD_OUT, chi);
  Bond bdo_D(BD_OUT, D);

  vector<Bond> bondsC;
  bondsC.push_back(bdi_chi);
  bondsC.push_back(bdo_chi);

  vector<Bond> bondsT;
  bondsT.push_back(bdi_D);
  bondsT.push_back(bdi_chi);
  bondsT.push_back(bdi_chi);

  vector<Bond> bondsAB;
  bondsAB.push_back(bdi_D);
  bondsAB.push_back(bdi_D);
  bondsAB.push_back(bdi_D);
  bondsAB.push_back(bdi_D);

  UniTensor<uni10_double64> C0(bondsC);
  UniTensor<uni10_double64> C1(bondsC);
  UniTensor<uni10_double64> C2(bondsC);
  UniTensor<uni10_double64> C3(bondsC);

  C0.randomize('U', 0, 1, 777);
  C0 *= 0.5;
  C1 = C0;
  C2 = C0;
  C3 = C0;

  UniTensor<uni10_complex128> C0C(C0);

  UniTensor<uni10_double64> Ta0(bondsT);
  UniTensor<uni10_double64> Ta1(bondsT);
  UniTensor<uni10_double64> Ta2(bondsT);
  UniTensor<uni10_double64> Ta3(bondsT);

  Ta0.randomize('U', 0, 1, 777);
  Ta0 *= 0.5;
  Ta1 = Ta0;
  Ta2 = Ta0;
  Ta3 = Ta0;

  UniTensor<uni10_complex128> Ta0C(Ta0);

  UniTensor<uni10_double64> Tb0(bondsT);
  UniTensor<uni10_double64> Tb1(bondsT);
  UniTensor<uni10_double64> Tb2(bondsT);
  UniTensor<uni10_double64> Tb3(bondsT);

  Tb0.randomize('U', 0, 1, 777);
  Tb0 *= 0.5;
  Tb1 = Tb0;
  Tb2 = Tb0;
  Tb3 = Tb0;

  UniTensor<uni10_complex128> Tb0C(Tb0);

  UniTensor<uni10_double64> A1(bondsAB);
  UniTensor<uni10_double64> A2(bondsAB);
  UniTensor<uni10_double64> B1(bondsAB);
  UniTensor<uni10_double64> B2(bondsAB);

  A1.randomize('U', 0, 1, 777);
  B1.randomize('U', 0, 1, 777);
  A1 *= 0.5;
  B1 *= 0.5;
  A2 = A1;
  B2 = B1;

  Network_dev cmt_net("ctm.net");
  C0C.setName("C0C");
  Ta0C.setName("Ta0C");

  UniTensor<uni10_double64> Tout;
  contract_args(Tout, cmt_net, C0, C1, C2, C3, Ta0, Ta1, Ta2, Ta3, Tb0, Tb1, Tb2, Tb3, A1, A2, B1, B2);
  cout << Tout;

  Network<uni10_double64> cmt_net_old("ctm_old.net");

  UniTensor<uni10_double64> Tout_old;
  contract_args(Tout_old, cmt_net_old, C0, C1, C2, C3, Ta0, Ta1, Ta2, Ta3, Tb0, Tb1, Tb2, Tb3, A1, A2, B1, B2);
  cout << Tout_old;

  C0.randomize();
  C0C = C0;

  UniTensor<uni10_complex128> ToutC;

  contract_args(ToutC, cmt_net, C0C, C1, C2, C3, Ta0, Ta1, Ta2, Ta3, Tb0C, Tb0C, Tb2, Tb3, A1, A2, B1, B2);
  cout << cmt_net;
  cout << ToutC;

  contract_args(Tout_old, cmt_net_old, C0, C1, C2, C3, Ta0, Ta1, Ta2, Ta3, Tb0, Tb1, Tb2, Tb3, A1, A2, B1, B2);

  cout << cmt_net_old.launch();

  return 0;

}
