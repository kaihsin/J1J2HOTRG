#include "5pess_tools/5pess_kagome.h"
#include "../../../hamiltonian/hamiltonian.h"

using namespace std;
using namespace uni10;

// It is a simple example for calculating the groud state energe of Ising model by itebd in 1 dimensional system.
// The tensors utilized in this example are without considering any symmetry.
//  
int main(){

  uni10_create();
  uni10_print_env_info();

  pess5_paras paras;
  paras.load_pess5_paras();

  UniTensor<uni10_double64> hamiltonian_d;
  UniTensor<uni10_complex128> hamiltonian_c;

  bool is_real = load_hamiltonian(hamiltonian_d, hamiltonian_c);

  Network_dev theta_net("5pess_net/theta.net");
  Network_dev state_net("5pess_net/state.net");
  Network_dev measure_net("5pess_net/measure.net");

  map< string, Network_dev* > net_list;
  net_list["theta"] = &theta_net;
  net_list["state"] = &state_net;
  net_list["measure"] = &measure_net;

  if(is_real){
    PESS_5<uni10_double64> pess_5_run(hamiltonian_d, paras, net_list);
    pess_5_run.Optimize();
  }
  else{
    PESS_5<uni10_complex128> pess_5_run(hamiltonian_c, paras, net_list);
    pess_5_run.Optimize();
  }

  uni10_destroy();

  return 0;
}
