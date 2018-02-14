#include "itebd_tools.h"
#include "itebd_1d.h"

//~/GitRepo/tensorlib/uni10/example/non-symmetry/return_linalg_funcs/itebd_1d/itebd_tools
template<typename T>
iTEBD_1D<T>::iTEBD_1D(const UniTensor<T>& _H,  const itebd_paras& paras): H(_H){

  fprintf(stdout, "\n");
  paras.print_info();
  fprintf(stdout, "\n");

  std::cout << H;

  dim = H.bond(0).dim();
  D = paras.D;
  max_N   = paras.max_N;
  tau     = paras.tau;
  eps     = paras.eps;
  cut_off = paras.cut_off;
  measure_per_n_iter = paras.measure_per_n_iter;

  max_D = (uni10_int)cut_off == -1 ? D : paras.max_D;
  this->init();

}

template<typename T>
void iTEBD_1D<T>::init(){

  vector<Bond> gamma_bds(3);
  gamma_bds[0] = Bond( BD_IN, D);
  gamma_bds[1] = Bond( BD_IN, D);
  gamma_bds[2] = Bond( BD_OUT, dim);
  gammas = vector< UniTensor<T> > ( 2, UniTensor<T>(gamma_bds));
  lambdas = vector<Matrix<T> > ( 2, Matrix<T>( D, D, true));

  ///create random matrix
  Matrix<T> temp = gammas[0].getBlock();
  uni10_rand( temp, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  Matrix<T> temp1(D, D, true);
  uni10_rand( temp1, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  for (uni10_int32 i=0; i!=2; i++){
    gammas[i].putBlock(temp);
    lambdas[i] = temp1;
  }

}

template<typename T>
iTEBD_1D<T>::~iTEBD_1D(){


}

// Set hamiltonian in iTEBD algorithm.
template<typename T>
void iTEBD_1D<T>::setHamiltonian(const UniTensor<T>& _H){

  H = _H;

}

// Get the gate 
template<typename T>
UniTensor<T> iTEBD_1D<T>::get_gate(const UniTensor<T>& _H){

  UniTensor<T> gate( _H.bond() );
  gate.putBlock( exph( -1.0*tau, H.const_getBlock()) );
  return gate;

}

template<typename T>
void iTEBD_1D<T>::Optimize(){

  fprintf(stdout, "Updating the matrix product states: \n\n");
  progressbar(0, 0, max_N, true);

  uni10_int32 ham_label[] = {4, 5, 6, 7};
  UniTensor<T> gate = get_gate(H);
  H.setLabel(ham_label);
  gate.setLabel(ham_label);

  for (uni10_uint64 i=0; i < max_N; i++){
    uni10_int32 i_l = i%2, i_r = (i+1)%2;

    ///contract theta
    bondcat( gammas[i_l], lambdas[i_r], 0);
    bondcat( gammas[i_l], lambdas[i_l], 1);
    bondcat( gammas[i_r], lambdas[i_r], 1);

    uni10_int32 gl_label[] = {1, 2, 4};
    uni10_int32 gr_label[] = {2, 3, 5};

    gammas[i_l].setLabel(gl_label);
    gammas[i_r].setLabel(gr_label);
    
    uni10_int32 theta_label[] = {1, 6, 3, 7};
    UniTensor<T> theta; 
    theta = contract( gammas[i_l], gammas[i_r], false );
    theta = contract( theta, gate, false ); ///now theta label 1, 3; 6, 7
    theta = permute( theta, theta_label, 2 );

    ///svd and update
    vector<Matrix<T> > usv = svd(theta.getBlock());
    resize( usv[0], usv.at(0).row(),            D, INPLACE);
    resize( usv[1],           D    ,            D, INPLACE);
    resize( usv[2],           D    , usv[2].col(), INPLACE);
    lambdas[i_l] = usv[1];
    lambdas[i_l] *= 1.0/norm( lambdas[i_l]);

    uni10_int32 gl1_label[] = {1, 4, 2};
    uni10_int32 gr1_label[] = {2, 3, 5};
    gammas[i_l] = permute(gammas[i_l], gl1_label, 2); 
    gammas[i_r] = permute(gammas[i_r], gr1_label, 1);

    gammas[i_l].putBlock( usv[0]); 
    gammas[i_r].putBlock( usv[2]);

    uni10_int32 gl2_label[] = {1, 2, 4};
    uni10_int32 gr2_label[] = {2, 3, 5};
    gammas[i_l] = permute( gammas.at(i_l), gl2_label, 2); 
    gammas[i_r] = permute( gammas.at(i_r), gr2_label, 2);

    bondrm( gammas[i_l], lambdas[i_r], 0);
    bondrm( gammas[i_r], lambdas[i_r], 1);

    ///measure
    uni10_int32 thetaH_label[] = {1, 3, 4, 5};
    if (i%measure_per_n_iter==0){
      vector<UniTensor<T> > gamma_now = gammas;
      bondcat( gamma_now[i_l], lambdas[i_r], 0);
      bondcat( gamma_now[i_l], lambdas[i_l], 1);
      bondcat( gamma_now[i_r], lambdas[i_r], 1);
      theta = contract( gamma_now[0], gamma_now[1], false); ///now theta label 1, 4; 3, 5
      UniTensor<T> theta_T =  theta;
      UniTensor<T> theta_H = contract( theta, H, false ); ///now theta_H label 1, 3; 6, 7
      theta_H.setLabel( thetaH_label );
      T norm  = contract( theta, theta_T, false )[0];
      T expec = contract( theta, theta_H, false )[0];
      progressbar(i+measure_per_n_iter, 0, max_N);
      cout.precision(8);
      cout.setf(ios::fixed, ios::floatfield);
      cout <<  ", ge: " << expec/norm  << "\r";
      std::cout.flush();

    }
    
  }

  fprintf(stdout, "\n\n");

}

template class iTEBD_1D<uni10_double64>;
template class iTEBD_1D<uni10_complex128>;
