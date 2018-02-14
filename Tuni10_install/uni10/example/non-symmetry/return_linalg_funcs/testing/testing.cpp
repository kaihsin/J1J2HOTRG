#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;


int main(){

  uni10_create();
  uni10_print_env_info();

  Matrix<uni10_double64> a(6, 6);
  uni10_rand( a, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  Matrix<uni10_complex128> b(6, 6);
  uni10_rand( b, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  Matrix<uni10_double64> ad(6, 6, true);
  uni10_rand( ad, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  Matrix<uni10_complex128> bd(6, 6, true);
  uni10_rand( bd, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);

  Matrix<uni10_double64> a2(a);

  vector< pair<void*, int> > pa(1, pair<void*, int>( (void*)&a, 1));
  pair<void*, int> pb = pa[0];

  cout << pa[0].first << "  " << pa[0].second << endl;
  cout << pb.first << "  " << pb.second << endl;

  cout << pa.size() << endl;
  pa[0] = pair<void*, int>((void*)&b, 2);
  cout << pa.size() << endl;

  cout << pa[0].first << "  " << pa[0].second << endl;
  cout << pb.first << "  " << pb.second << endl;


  vector< pair<void*, int> > mlist(4);
  mlist[0] = pair<void*, int>((void*)&b,  b.typeID());
  mlist[1] = pair<void*, int>((void*)&b,  b.typeID());
  mlist[2] = pair<void*, int>((void*)&bd, bd.typeID());
  mlist[3] = pair<void*, int>((void*)&bd, bd.typeID());

  vector<Matrix<double>*> mlist_d(4);
  mlist_d[0] = &a;
  mlist_d[1] = &a;
  mlist_d[2] = &ad;
  mlist_d[3] = &ad;

  cout << dot(a, a);
  cout << dot(bd, bd);
  cout << dot(dot(dot(dot(dot(dot(a, bd),dot(b, ad)), a), a), b), b);
  //cout << dot( dot( dot( dot(a, a), dot(ad, ad) ), a), a );
  //cout << dot(dot(dot(b, b), dot(bd, bd)), dot( dot(a, ad), a ));

  Matrix<uni10_complex128> mout_z_args;
  Matrix<uni10_double64> mout_d_args;
  //dot_args(mout_d_args, a, a, ad, ad, a, a);
  dot_args(mout_z_args, a, bd, b, ad, a, a, b, b);
  cout << mout_z_args;

  return 0;

}
