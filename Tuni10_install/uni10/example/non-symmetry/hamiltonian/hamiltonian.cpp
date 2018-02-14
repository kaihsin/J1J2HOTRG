#include "hamiltonian.h"

uni10::Bond spin_bond(float spin, uni10::bondType btype, bool U1){

  spin_check(spin);
  int dim = spin * 2 + 1;

  if(U1){
    std::vector<uni10::Qnum> qnums(dim);
    int halfint = true;
    if(spin == floor(spin))
      halfint = false;
    for(int i = 0; i < dim; i++){
      int s = spin - i;
      if(halfint){
        s = spin + 0.5 - i;
        if(s <= 0)
          s--;
      }
      qnums[i] = uni10::Qnum(s);
    }
    return uni10::Bond(btype, qnums);
  }
  else{

    return uni10::Bond(btype, dim);

  }

}

uni10::UniTensor<double> XXZ(float Jx, float Jz, float spin){

  uni10::Matrix<double> sp = matSp(spin);
  uni10::Matrix<double> sm = matSm(spin);
  uni10::Matrix<double> sz = matSz(spin);
  uni10::Matrix<double> ham = (double)Jz*uni10::otimes(sz, sz);
  ham += Jx* 0.5 * (uni10::otimes(sp, sm) + uni10::otimes(sm, sp));
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  std::vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor<double> H(bonds, "Heisenberg");
  H.putBlock(ham);
  
  return H;

}

uni10::UniTensor<double> Heisenberg(float spin, double J){

  uni10::Matrix<double> sp = matSp(spin);
  uni10::Matrix<double> sm = matSm(spin);
  uni10::Matrix<double> sz = matSz(spin);
  uni10::Matrix<double> ham = uni10::otimes(sz, sz);
  ham += 0.5 * (uni10::otimes(sp, sm) + uni10::otimes(sm, sp));
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  std::vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor<double> H(bonds, "Heisenberg");
  H.putBlock(ham);
  return J * H;

}

uni10::UniTensor<double> Heisenberg_U1(float spin, double J){

}

uni10::UniTensor<double> transverseIsing(float spin, float h, bool isAnti){

  uni10::Matrix<double> sx = matSx(spin);
  uni10::Matrix<double> sz = matSz(spin);
  uni10::Matrix<double> I(sx.row(), sx.col(), true);
  I.identity();
  uni10::Matrix<double> ham = (isAnti) ? uni10::otimes((double)2*sz, (double)2*sz) : ((double)-1.) * uni10::otimes((double)2*sz, (double)2*sz); // otimes(sigma_z, sizga_z);
  uni10::Matrix<double> sxl = uni10::otimes((h/(double)2) * (double)2*sx, I);
  uni10::Matrix<double> sxr = uni10::otimes(I, (h/(double)2) * (double)2*sx);
  ham = ham + 0.5 * (sxl + sxr) ;
  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  std::vector<uni10::Bond> bonds(2, bdi);
  bonds.push_back(bdo);
  bonds.push_back(bdo);
  uni10::UniTensor<double> H(bonds, "transverseIsing");
  H.putBlock(ham);
  return H;

}

uni10::UniTensor<double> theModel(float spin, int i, double delta, double Dz, double hz, double dx){

}

uni10::UniTensor<double> JQmodel(double J, double Q){

}

uni10::UniTensor<double> periodicHamiltonian(int N, const uni10::UniTensor<double>& H0){

}

uni10::UniTensor<double> J1J2model(double J1, double J2){

}

uni10::UniTensor<double> directSum(size_t s1, size_t s2, const uni10::UniTensor<double>& H2, const uni10::UniTensor<double>& H_all){

}

bool load_hamiltonian(uni10::UniTensor<double>& hamiltonian_d, uni10::UniTensor<std::complex<double> >& hamiltonian_c, std::string fname){

  bool is_real = false;
  char name[256];
  double para1, para2, para3, para4;
  

  FILE* rcfp = fopen(".hamrc", "r");
  int max_len = 256;
  char buffer[max_len];

  char* pch;
  while(fgets(buffer, max_len, rcfp)){

    buffer[strlen(buffer)-1] = '\0';

    pch = strtok(buffer, ": \n");
    pch = strtok(NULL, ": \n");

    if(strcmp ("HAMILTONIAN", buffer)==0)
      strcpy(name, pch);

    else if(strcmp ("PARA1", buffer)==0)
      para1 = atof(pch);

    else if(strcmp ("PARA2", buffer)==0)
      para2 = atof(pch);

    else if(strcmp ("PARA3", buffer)==0)
      para3 = atof(pch);

    else if(strcmp ("PARA4", buffer)==0)
      para4 = atof(pch);

    else if(buffer[0] =='#' || pch == NULL)
      continue;

    else{
      fprintf(stdout, "%s", "Setting the parameters with wrong names.");
      exit(0);
    }

  }

  if(std::string(name) == "transverseIsing"){
    is_real = true;
    hamiltonian_d = transverseIsing(para1, para2, para3);
  }

  else if(std::string(name) == "Heisenberg"){
    is_real = true;
    hamiltonian_d = Heisenberg(para1, para2);
  }else{

    uni10_error_msg(true, "%s", "Unexpected error. Please connect to the developers of uni10.");
  
  }

  fclose(rcfp);

  return is_real;

};


