#include "operator.h"

using namespace std;
using namespace uni10;

int main(){

  cout << "/****************************/" << endl;
  cout << "/*           Sx             */" << endl;
  cout << "/****************************/" << endl;
  cout << matSx(0.5);
  
  cout << "/****************************/" << endl;
  cout << "/*           Sz             */" << endl;
  cout << "/****************************/" << endl;
  cout << matSz(0.5);

  cout << "/****************************/" << endl;
  cout << "/*           Sp             */" << endl;
  cout << "/****************************/" << endl;
  cout << matSp(0.5);

  cout << "/****************************/" << endl;
  cout << "/*           Sm             */" << endl;
  cout << "/****************************/" << endl;
  cout << matSm(0.5);

  return 0;
}
