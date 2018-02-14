#include "uni10/uni10.hpp"

#include "iostream"
#include "string"

using namespace std;
using namespace uni10;


int main(){

  int chi = 5;
  int D = 4;

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

  UniTensor<uni10_complex128> C(bondsC);


  UniTensor<uni10_complex128> CC(bondsC);

  UniTensor<uni10_complex128> Ta(bondsT);
  UniTensor<uni10_complex128> Tb(bondsT);
  UniTensor<uni10_complex128> A(bondsAB);
  UniTensor<uni10_complex128> B(bondsAB);

  CC.randomize();
  C.randomize();
  Ta.randomize();
  Tb.randomize();
  A.randomize();
  B.randomize();

  C *= (1/2.);
  Ta *= (1/2.);
  Tb *= (1/2.);
  A *= (1/2.);
  B *= (1/2.);

  Network_dev cmt_net("ctm.net");

  cmt_net.putTensor(0, C);
  cmt_net.putTensor(1, C);
  cmt_net.putTensor(2, C);
  cmt_net.putTensor(3, C);

  cmt_net.putTensor(4, Ta);
  cmt_net.putTensor(5, Ta);
  cmt_net.putTensor(6, Ta);
  cmt_net.putTensor(7, Ta);

  cmt_net.putTensor(8, Tb);
  cmt_net.putTensor(9, Tb);
  cmt_net.putTensor(10, Tb);
  cmt_net.putTensor(11, Tb);

  cmt_net.putTensor(12, A);
  cmt_net.putTensor(13, A);
  cmt_net.putTensor(14, B);
  cmt_net.putTensor(15, B);

  UniTensor<uni10_complex128> Tout;

  cmt_net.pre_construct();
  cout << cmt_net;
  cout << cmt_net.get_contract_order() << endl;
  
  cmt_net.launch(Tout);
  cout << Tout;

  Network<uni10_complex128> cmt_net_old("ctm_old.net");

  cmt_net_old.putTensor(0, C);
  cmt_net_old.putTensor(1, C);
  cmt_net_old.putTensor(2, C);
  cmt_net_old.putTensor(3, C);

  cmt_net_old.putTensor(4, Ta);
  cmt_net_old.putTensor(5, Ta);
  cmt_net_old.putTensor(6, Ta);
  cmt_net_old.putTensor(7, Ta);

  cmt_net_old.putTensor(8, Tb);
  cmt_net_old.putTensor(9, Tb);
  cmt_net_old.putTensor(10, Tb);
  cmt_net_old.putTensor(11, Tb);

  cmt_net_old.putTensor(12, A);
  cmt_net_old.putTensor(13, A);
  cmt_net_old.putTensor(14, B);
  cmt_net_old.putTensor(15, B);

  cout << cmt_net_old.launch();

  /*
  C.randomize();
  Ta.randomize();
  Tb.randomize();
  A.randomize();
  B.randomize();

  C *= (1/2.);
  Ta *= (1/2.);
  Tb *= (1/2.);
  A *= (1/2.);
  B *= (1/2.);

  Network_dev cmt_net_odd("ctm_odd.net");
  cmt_net_odd.putTensor(0, C);
  cmt_net_odd.putTensor(1, C);
  cmt_net_odd.putTensor(2, C);
  cmt_net_odd.putTensor(3, C);

  cmt_net_odd.putTensor(4, Ta);
  cmt_net_odd.putTensor(5, Ta);
  cmt_net_odd.putTensor(6, Ta);
  cmt_net_odd.putTensor(7, Ta);

  cmt_net_odd.putTensor(8, Tb);
  cmt_net_odd.putTensor(9, Tb);
  cmt_net_odd.putTensor(10, Tb);
  cmt_net_odd.putTensor(11, Tb);

  cmt_net_odd.putTensor(12, A);
  cmt_net_odd.putTensor(13, A);
  cmt_net_odd.putTensor(14, B);

  cmt_net_odd.launch(Tout, "QQ");
  cout << cmt_net_odd;
  cout << Tout;


  Network<uni10_complex128> cmt_net_old_odd("ctm_old_odd.net");
  cmt_net_old_odd.putTensor(0, C);
  cmt_net_old_odd.putTensor(1, C);
  cmt_net_old_odd.putTensor(2, C);
  cmt_net_old_odd.putTensor(3, C);

  cmt_net_old_odd.putTensor(4, Ta);
  cmt_net_old_odd.putTensor(5, Ta);
  cmt_net_old_odd.putTensor(6, Ta);
  cmt_net_old_odd.putTensor(7, Ta);

  cmt_net_old_odd.putTensor(8, Tb);
  cmt_net_old_odd.putTensor(9, Tb);
  cmt_net_old_odd.putTensor(10, Tb);
  cmt_net_old_odd.putTensor(11, Tb);

  cmt_net_old_odd.putTensor(12, A);
  cmt_net_old_odd.putTensor(13, A);
  cmt_net_old_odd.putTensor(14, B);

  cout << cmt_net_old_odd.launch("QQ");
  */

  return 0;

}
