#ifndef _MPO_HAMILTONIAN_H_
#define _MPO_HAMILTONIAN_H_

#include "../operator/operator.h"
#include "../hamiltonian/hamiltonian.h"

//======================================

class MPO {

public:
	/// constructor
	MPO( int dim, char loc );

	void putTensor( uni10::UniTensor<double> op, int row_idx, int col_idx );
	uni10::UniTensor<double> launch();

private:
	int virt_dim;
	char mpo_loc;
	uni10::UniTensor<double> mpo_frame;
	uni10::UniTensor<double> mpo;

};

//======================================

std::vector<uni10::UniTensor<double>> mpoXXZ(float Jx, float Jz, float spin = 0.5);
std::vector<uni10::UniTensor<double>> mpoITF(float h, float spin = 0.5);

bool load_ham_mpo(std::vector<uni10::UniTensor<double>>& ham_mpo_d, std::string fname = "");

#endif
