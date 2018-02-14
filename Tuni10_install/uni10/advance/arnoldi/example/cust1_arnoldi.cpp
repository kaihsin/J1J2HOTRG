#include "uni10/uni10.hpp"

#include "../../../example/nsy_mpo/mpo.h"

#include "../arnoldi.h"

using namespace std;
using namespace uni10;

int main(){


  UniTensor<uni10_double64> WL = TranLongIsing_MPO(1., 1., 1.);
  UniTensor<uni10_double64> WR = WL;
  UniTensor<uni10_double64> EnvL = Edge_MPO(3, 2, "L");
  UniTensor<uni10_double64> EnvR = Edge_MPO(3, 0, "R");
  
  vector<Bond> bondsGS;
  bondsGS.push_back(Bond(BD_IN, 1));
  bondsGS.push_back(Bond(BD_IN, 2));
  bondsGS.push_back(Bond(BD_IN, 2));
  bondsGS.push_back(Bond(BD_IN, 1));
  UniTensor<uni10_double64> GS(bondsGS);

  Network<uni10_double64> l_net("../../../example/nsy_dmrg_1d/Diagrams/Lanczos.net");

  Matrix<uni10_double64> V(4, 1);
  uni10_rand(V, uni10_mt19937, uni10_uniform_real, -1, 1, uni10_clock);
  double Norm = norm( V );
  V *= (1./Norm);
  GS.putBlock(V);

  Network<uni10_double64> test("test.net");
  test.putTensor("EnvL", EnvL);
  test.putTensor("WL", WL);
  test.putTensor("WR", WR);
  test.putTensor("EnvR", EnvR);
  UniTensor<uni10_double64> tmp = test.launch();
  vector<Matrix<uni10_double64> > EigH = eigh(tmp.getBlock());

  cout << EigH[0] << endl;
  cout << EigH[1] << endl;

  uni10_arnoldi_custom1_paras<uni10_double64> cust1_paras(EnvL, WL, WR, EnvR, l_net, V);
  uni10_arnoldi_paras<uni10_double64, uni10_arnoldi_custom1_paras<uni10_double64> > l_paras(cust1_paras);

  uni10_uint64 iter=0;

  int info;
  uni10_complex128 E = Arnoldi(l_paras, iter, 10, 2, info);

  if(info != -1)
    cout << "  Not converge " << endl;

  cout << "E: " << E << endl;
  GS.putBlock(V);
  cout << GS;

  return 0;

}
