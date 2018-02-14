#include "../../mpo_hamiltonian/mpo_hamiltonian.h"
#include "idmrg_tools/idmrg_1d.h"
#include "idmrg_tools/mpsinf.h"

using namespace std;
using namespace uni10;

// It is a simple example for calculating the groud state energe of Ising model by idmrg in 1 dimensional system.
// The tensors utilized in this example are without considering any symmetry.
//
int main(){

  uni10_create();
  uni10_print_env_info();

  idmrg_paras paras;
  paras.load_idmrg_paras();
  paras.print_info();

  std::vector<UniTensor<double> > ham_mpo_d;

  bool is_real = load_ham_mpo(ham_mpo_d);

  if(is_real){
    int d = ham_mpo_d[1].bond()[1].dim();
    MPSInf mps(paras.chi, d, 2);
    mps.randomize();
    mps.idmrg(ham_mpo_d, paras);
  }
  else{
    std::cerr << "Complex type not supported." << '\n';
  }

  uni10_destroy();

  return 0;
}
