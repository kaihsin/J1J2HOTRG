#include "itebd_tools.h"
#include "itebd_2d.h"

//~/GitRepo/tensorlib/uni10/example/non-symmetry/return_linalg_funcs/itebd_1d/itebd_tools
template<typename T>
iTEBD_2D<T>::iTEBD_2D(const UniTensor<T>& _H,  const itebd_paras& paras): H(_H){

  fprintf(stdout, "\n");
  paras.print_info();
  fprintf(stdout, "\n");

  cout << H;

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
void iTEBD_2D<T>::init(){

  vector<Bond> gamma_bds( 4, Bond( BD_OUT, D) );
  gamma_bds.insert( gamma_bds.begin(), Bond( BD_IN, 2 ) );
  gammas = vector< UniTensor<T> > ( 2, UniTensor<T>(gamma_bds));
  lambdas = vector<Matrix<T> > ( 4, Matrix<T>( D, D, true));

  ///create random matrix
  Matrix<T> temp = gammas[0].getBlock();
  uni10_rand( temp, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  Matrix<T> temp1(D, D, true);
  uni10_rand( temp1, uni10_mt19937, uni10_normal, 0, 1, uni10_clock);
  for (uni10_int32 i=0; i!=gammas.size(); i++){
    gammas[i].putBlock(temp);
  }
  for (uni10_int32 i=0; i!=lambdas.size(); i++){
    lambdas[i] = temp1;
  }

}

template<typename T>
iTEBD_2D<T>::~iTEBD_2D(){


}

// Set hamiltonian in iTEBD algorithm.
template<typename T>
void iTEBD_2D<T>::setHamiltonian(const UniTensor<T>& _H){

  H = _H;

}

// Get the gate 
template<typename T>
UniTensor<T> iTEBD_2D<T>::get_gate(const UniTensor<T>& _H){

  UniTensor<T> gate( _H.bond() );
  gate.putBlock( exph( -1.0*tau, H.const_getBlock()) );
  return gate;

}

template<typename T>
void iTEBD_2D<T>::Optimize(){

  fprintf(stdout, "Updating the matrix product states: \n\n");
  progressbar(0, 0, max_N, true);

  UniTensor<T> gate = get_gate(H);

  for (uni10_uint64 i=0; i < max_N; i++){
    for ( int idir=0; idir!=4; idir++){
      bondcatAll( idir );
      UniTensor<T> theta = contractGammas( idir );
      contractTwoSiteOp( theta, gate );

      ///svd and update
      vector<Matrix<T> > usv = svd(theta.getBlock());
      resize( usv[0], usv.at(0).row(),            D, INPLACE);
      resize( usv[1],               D,            D, INPLACE);
      resize( usv[2],               D, usv[2].col(), INPLACE);

      //update gammas and lambda
      lambdas[idir]  = usv[1];
      lambdas[idir] *= 1.0/norm( lambdas[idir] );

      updateGammas( usv[0], usv[2], idir );

      //recover
      bondrmAll( idir );
      bondcat( gammas[0], lambdas[idir], (idir+2)%4+1 );
    }

    //measure
    if (i%measure_per_n_iter==0){
      progressbar(i+measure_per_n_iter, 0, max_N);
      T average = 0;
      for ( int idir=0; idir!=4; idir++ ){
        T norm = measureNorm( idir );
        T expectation = measureExpe( H, idir );
        //printf( "%12.4e%2s", expectation/norm, "" );
        average += expectation/norm;
      }
      average *= 0.25;
      cout.precision(8);
      cout.setf(ios::fixed, ios::floatfield);
      cout <<  ", ge: " << average  << "\r";
      std::cout.flush();
    }
  }

  fprintf(stdout, "\n\n");

}

template<typename T>
void iTEBD_2D<T>::bondcatAll( const int idir){
  uni10_int32 i_l =0, i_r =1;
  enum{ bond_phy, bond_l, bond_u, bond_r, bond_d};//left,up,right,down    

  for( int i=bond_l; i<= bond_d; i++){
    bondcat( gammas[i_l], lambdas[i-1], i);
  }
  for( int j=bond_l; j<= bond_d; j++){
    if ( (idir+1)!=j){
      bondcat( gammas[i_r], lambdas[j-1], j);
    }
    else {}
  }
}

template<typename T>
void iTEBD_2D<T>::bondrmAll( const int idir ){
  uni10_int32 i_l =0, i_r =1;
  enum{ bond_phy, bond_l, bond_u, bond_r, bond_d};//left,up,right,down    

  for( int i=bond_l; i<= bond_d; i++){
    bondrm( gammas[i_l], lambdas[i-1], i);
  }
  for( int j=bond_l; j<= bond_d; j++){
    if ( (idir+1)!=j){
      bondrm( gammas[i_r], lambdas[j-1], j);
    }
    else {}
  }
}

/*
template<typename T>
void iTEBD_2D<T>::bondrmSix( const int idir ){
  uni10_int32 i_l =0, i_r =1;
  enum{ bond_phy, bond_l, bond_u, bond_r, bond_d};//left,up,right,down    

  for( int i=bond_l; i<= bond_d; i++){
    bondrm( gammas[i_l], lambdas[i-1], i);
  }
  for( int j=bond_l; j<= bond_d; j++){
    if ((idir+1)!=j){
      bondrm( gammas[i_r], lambdas[j-1], j);
    }
    else {}
  }
}
*/

template<typename T>
UniTensor<T> iTEBD_2D<T>::contractGammas( const int idir ){
  vector<int> gamma0Lab = { 1, 2, 3, 4, 5};
  vector<int> gamma1Lab = { 6, 7, 8, 9, 10};
  gamma0Lab.at( (idir+2)%4+1 ) = 0;
  gamma1Lab.at(  idir+1 ) = 0;
  gammas[0].setLabel( gamma0Lab );
  gammas[1].setLabel( gamma1Lab );
  UniTensor<T> theta = contract( gammas[0], gammas[1], false );
  return theta;
}

template<typename T>
void iTEBD_2D<T>::contractTwoSiteOp( UniTensor<T> &theta, UniTensor<T> &twoSiteOp ){
  vector<int> oldLab = theta.label();
  const int inbdn = theta.inBondNum();
  twoSiteOp.setLabel( {1, 6, -1, -6} );
  theta = contract( theta, twoSiteOp, false );
  vector<int> newLab = theta.label();
  newLab.at(6) = 1;
  newLab.at(7) = 6;
  theta.setLabel( newLab );
  theta = permute( theta, oldLab, inbdn );
}

template<typename T>
void iTEBD_2D<T>::updateGammas( Matrix<T> &u, Matrix<T> &vT, const int idir ){
  vector<Bond> newLbds( 3, Bond( BD_IN, D ) );
  newLbds.insert( newLbds.begin(), Bond( BD_IN, dim ) );
  newLbds.push_back( Bond( BD_OUT, D ) );
  vector<Bond> newRbds( 3, Bond( BD_OUT, D ) );
  newRbds.insert( newRbds.begin(), Bond( BD_OUT, dim) );
  newRbds.insert( newRbds.begin(), Bond( BD_IN, D ) );
  UniTensor<T> newL( newLbds );
  UniTensor<T> newR( newRbds );
  newL.putBlock( u );
  newR.putBlock( vT );

  vector<int> newLlabs = { -1, 1, 2, 3, 0 };
  newL.setLabel( newLlabs );
  int posit = (idir+2)%4+1;
  newLlabs.insert( newLlabs.begin()+posit, newLlabs.back() );
  newLlabs.pop_back();
  gammas[0] = permute( newL, newLlabs, 1 );

  vector<int> newRlabs = { 0, -1, 1, 2, 3 };
  newR.setLabel( newRlabs );
  newRlabs.insert( newRlabs.begin()+2+idir, newRlabs.at(0) );
  newRlabs.erase( newRlabs.begin() );
  gammas[1] = permute( newR, newRlabs, 1 );
}

template<typename T>
T iTEBD_2D<T>::measureNorm( const int idir ){
  bondcatAll( idir );
  UniTensor<T> theta = contractGammas( idir );
  UniTensor<T> thetaT = theta;
  thetaT = transpose( dagger( theta ) );//
  UniTensor<T> out = contract( theta, thetaT, false );
  bondrmAll( idir );
  return out[0];
}

template<typename T>
T iTEBD_2D<T>::measureExpe( UniTensor<T> &twoSiteOp, const int idir ){
  bondcatAll( idir );
  UniTensor<T> theta = contractGammas( idir );
  UniTensor<T> thetaT = theta;
  thetaT = transpose( dagger( theta ) );//
  contractTwoSiteOp( theta, twoSiteOp );
  UniTensor<T> out = contract( theta, thetaT, false );
  bondrmAll( idir );
  return out[0];
}

template class iTEBD_2D<uni10_double64>;
template class iTEBD_2D<uni10_complex128>;
