#include "../../hamiltonian/hamiltonian.h"
#include "itebd_tools/itebd_1d.h"

using namespace std;
using namespace uni10;

// It is a simple example for calculating the groud state energe of Ising model by itebd in 1 dimensional system.
// The tensors utilized in this example are without considering any symmetry.
//  
int main(){

  uni10_create();
  uni10_print_env_info();

  itebd_paras paras;
  paras.load_itebd_paras();

  UniTensor<uni10_double64> hamiltonian_d;
  UniTensor<uni10_complex128> hamiltonian_c;

  bool is_real = load_hamiltonian(hamiltonian_d, hamiltonian_c);

  if(is_real){
    iTEBD_1D<uni10_double64> itebd_run(hamiltonian_d, paras);
    itebd_run.Optimize();
  }
  else{
    iTEBD_1D<uni10_complex128> itebd_run(hamiltonian_c, paras);
    itebd_run.Optimize();
  }

  uni10_destroy();

  return 0;
}
