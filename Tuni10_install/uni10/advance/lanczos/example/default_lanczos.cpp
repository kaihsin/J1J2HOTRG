#include "uni10/uni10.hpp"

#include "../../../example/nsy_mpo/mpo.h"

#include "../customized_structure.h"
#include "../lanczos.h"

using namespace std;
using namespace uni10;

int main(){

  uni10_uint64 Rnum = 20;
  uni10_uint64 Cnum = 20;

  Matrix<uni10_double64> M(Rnum, Cnum);
  uni10_rand(M, uni10_mt19937, uni10_normal, -1, 1, uni10_clock);

  // Symmetrize
  for(int i = 0; i < (int)M.row(); i++)
    for(int j = i+1; j < (int)M.col(); j++)
      M[j*M.col()+i] = M[i*M.col()+j];

  vector<Matrix<uni10_double64> > EigH = eigh(M);

  cout << "========= lapack ==========";

  cout << EigH[0] << endl;
  cout << EigH[1] << endl;

  cout << dot(transpose(EigH[1]), EigH[1]);


  cout << endl << endl;
  cout << "========= Lanczos1 ==========";
  cout << endl << endl;

  Matrix<uni10_double64> V(M.row(), 1);
  uni10_rand(V, uni10_mt19937, uni10_normal, 0, 1, 777);
  uni10_double64 Norm = norm( V );
  V *= (1./Norm);

  uni10_lanczos_default_paras<uni10_double64> default_paras(M, V);
  uni10_lanczos_paras<uni10_double64, uni10_lanczos_default_paras<uni10_double64> > l_paras(default_paras);

  uni10_uint64 iter=0;

  Matrix<uni10_double64> eigv;
   
  uni10_int32 info;
  Lanczos(eigv, 3, l_paras, iter, M.row(), 10, info);

  if(info != -1)
    cout << "  Not converge " << endl;

  cout << eigv << endl;
  cout << V;

  cout << dot(transpose(V), V);

  cout << endl << endl;
  cout << "========= Lanczos2 ==========";
  cout << endl << endl;

  V.assign(M.row(), 1);
  uni10_rand(V, uni10_mt19937, uni10_normal, 0, 1, 777);
  Norm = norm( V );
  V *= (1./Norm);

  iter=0;

  uni10_lanczos_default_paras<uni10_double64> default2_paras(M, V);
  uni10_lanczos_paras<uni10_double64, uni10_lanczos_default_paras<uni10_double64> > l2_paras(default2_paras);
  uni10_double64 E = Lanczos(l2_paras, iter, M.row(), 10, info);

  cout << "E: " << E << endl;
  cout << V;

  return 0;

}
