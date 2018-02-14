#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_

#include "../operator/operator.h"

uni10::Bond spin_bond(float spin, uni10::bondType btype, bool U1=false);

uni10::UniTensor<double> XXZ(float Jx, float Jz, float spin = 0.5);

uni10::UniTensor<double> Heisenberg(float spin=0.5, double J=1.0);

uni10::UniTensor<double> Heisenberg_U1(float spin=0.5, double J=1.0);

uni10::UniTensor<double> transverseIsing(float spin, float h, bool isAnti);

uni10::UniTensor<double> theModel(float spin, int i, double delta, double Dz, double hz, double dx);

uni10::UniTensor<double> JQmodel(double J, double Q);

uni10::UniTensor<double> periodicHamiltonian(int N, const uni10::UniTensor<double>& H0);

uni10::UniTensor<double> J1J2model(double J1, double J2);

uni10::UniTensor<double> directSum(size_t s1, size_t s2, const uni10::UniTensor<double>& H2, const uni10::UniTensor<double>& H_all);

bool load_hamiltonian(uni10::UniTensor<double>& hamiltonian_d, uni10::UniTensor< std::complex<double> >& hamiltonian_c, std::string fname = "");

#endif
