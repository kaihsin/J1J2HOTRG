#include <string>

#include "mpo_hamiltonian.h"


//======================================

MPO::MPO( int dim, char loc ) {
  ///
  virt_dim = dim;
  mpo_loc = loc;

  if ( loc == 'l' ) {
    std::vector< uni10::Bond > bonds_l;
    bonds_l.push_back( uni10::Bond( uni10::BD_OUT, virt_dim ) );
    mpo_frame.assign( bonds_l );
  }
  else if ( loc == 'r' ) {
    std::vector< uni10::Bond > bonds_r;
    bonds_r.push_back( uni10::Bond( uni10::BD_IN, virt_dim ) );
    mpo_frame.assign( bonds_r );
  }
  else {
    std::vector< uni10::Bond > bonds_m;
    bonds_m.push_back( uni10::Bond( uni10::BD_IN, virt_dim ) );
    bonds_m.push_back( uni10::Bond( uni10::BD_OUT, virt_dim ) );
    mpo_frame.assign( bonds_m );
  }

  mpo_frame.set_zeros();
}

//======================================

void MPO::putTensor( uni10::UniTensor<double> op, int row_idx, int col_idx ) {
  ///
  uni10::UniTensor<double> location( mpo_frame.bond() );
  location.set_zeros();
  uni10::Matrix<double> temp = location.getBlock();
  temp[row_idx * temp.col() + col_idx] = 1.0;
  location.putBlock( temp );

  if ( mpo.bondNum() == 0 )
    mpo = uni10::otimes( location, op );
  else
    mpo = mpo + uni10::otimes( location, op );
}

//======================================

uni10::UniTensor<double> MPO::launch() {
  ///
  if ( mpo.bondNum() == 0 )
    return mpo_frame;
  else
    return mpo;
}

//======================================

std::vector<uni10::UniTensor<double>> mpoXXZ(float Jx, float Jz, float spin){
  ///
  std::vector< uni10::UniTensor<double> > ham;

  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  std::vector<uni10::Bond> bonds;
  bonds.push_back(bdi);
  bonds.push_back(bdo);
  uni10::UniTensor<double> Sp(bonds);
  uni10::UniTensor<double> Sm(bonds);
  uni10::UniTensor<double> Sz(bonds);
  uni10::UniTensor<double> Id(bonds);
  Sp.putBlock(matSp(spin));
  Sm.putBlock(matSm(spin));
  Sz.putBlock(matSz(spin));
  Id.identity();

  MPO mpo_l(5, 'l');
  MPO mpo_m(5, 'm');
  MPO mpo_r(5, 'r');

  mpo_l.putTensor( 0.5*Jx*Sm, 0, 1 );
  mpo_l.putTensor( 0.5*Jx*Sp, 0, 2 );
  mpo_l.putTensor(     Jz*Sz, 0, 3 );
  mpo_l.putTensor(        Id, 0, 4 );

  mpo_m.putTensor(        Id, 0, 0 );
  mpo_m.putTensor(        Sp, 1, 0 );
  mpo_m.putTensor(        Sm, 2, 0 );
  mpo_m.putTensor(        Sz, 3, 0 );
  mpo_m.putTensor( 0.5*Jx*Sm, 4, 1 );
  mpo_m.putTensor( 0.5*Jx*Sp, 4, 2 );
  mpo_m.putTensor(     Jz*Sz, 4, 3 );
  mpo_m.putTensor(        Id, 4, 4 );

  mpo_r.putTensor(        Id, 0, 0 );
  mpo_r.putTensor(        Sp, 1, 0 );
  mpo_r.putTensor(        Sm, 2, 0 );
  mpo_r.putTensor(        Sz, 3, 0 );

  ham.push_back(mpo_l.launch());
  ham.push_back(mpo_m.launch());
  ham.push_back(mpo_r.launch());

  return ham;
}

std::vector<uni10::UniTensor<double>> mpoITF(float h, float spin){
  ///
  std::vector< uni10::UniTensor<double> > ham;

  uni10::Bond bdi = spin_bond(spin, uni10::BD_IN);
  uni10::Bond bdo = spin_bond(spin, uni10::BD_OUT);
  std::vector<uni10::Bond> bonds;
  bonds.push_back(bdi);
  bonds.push_back(bdo);
  uni10::UniTensor<double> Sx(bonds);
  uni10::UniTensor<double> Sz(bonds);
  uni10::UniTensor<double> Id(bonds);
  Sx.putBlock(matSx(spin));
  Sz.putBlock(matSz(spin));
  Id.identity();

  MPO mpo_l(3, 'l');
  MPO mpo_m(3, 'm');
  MPO mpo_r(3, 'r');

  mpo_l.putTensor(   h*Sx, 0, 0 );
  mpo_l.putTensor( -1.*Sz, 0, 1 );
  mpo_l.putTensor(     Id, 0, 2 );

  mpo_m.putTensor(     Id, 0, 0 );
  mpo_m.putTensor(     Sz, 1, 0 );
  mpo_m.putTensor(   h*Sx, 2, 0 );
  mpo_m.putTensor( -1.*Sz, 2, 1 );
  mpo_m.putTensor(     Id, 2, 2 );

  mpo_r.putTensor(     Id, 0, 0 );
  mpo_r.putTensor(     Sz, 1, 0 );
  mpo_r.putTensor(   h*Sx, 2, 0 );

  ham.push_back(mpo_l.launch());
  ham.push_back(mpo_m.launch());
  ham.push_back(mpo_r.launch());

  return ham;
}

bool load_ham_mpo(std::vector<uni10::UniTensor<double>>& ham_mpo_d, std::string fname){

  bool is_real = false;
  char name[256];
  double para1, para2, para3, para4;

  FILE* rcfp = fopen(".hamrc", "r");
  int max_len = 256;
  char buffer[max_len];

  char* pch;
  while(fgets(buffer, max_len, rcfp)){

    pch = strtok(buffer, ":");
    pch = strtok(NULL, ":");

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

  std::string name_str = name;
  if( name_str.find("XXZ") != std::string::npos ){
    is_real = true;
    ham_mpo_d = mpoXXZ(para2, para3, para1);
  }
  if( name_str.find("ITF") != std::string::npos ){
    is_real = true;
    ham_mpo_d = mpoITF(para2, para1);
  }

  fclose(rcfp);

  return is_real;
};
