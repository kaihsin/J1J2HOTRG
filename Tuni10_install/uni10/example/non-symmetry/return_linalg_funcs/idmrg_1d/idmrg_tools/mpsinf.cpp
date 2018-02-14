#include "mpsinf.h"

/*========================================

  MPSInf Class functions

  ========================================*/

MPSInf::MPSInf(int X, int d, int L) {
  /// object constructor
  lat_size = L;	// lattice size
  dim_phys = d;	// physical dimension
  chi_max = X;	// maximum bond dimension
}

MPSInf::~MPSInf() {
  /// object destructor
}

uni10::UniTensor<double> MPSInf::initGamma(int chi1, int chi2, int d) {
  /// initialize a gamma tensor
  std::vector<uni10::Bond> bond_gam;
  bond_gam.push_back( uni10::Bond(uni10::BD_IN, chi1) );
  bond_gam.push_back( uni10::Bond(uni10::BD_IN, d) );
  bond_gam.push_back( uni10::Bond(uni10::BD_OUT, chi2) );
  return uni10::UniTensor<double>(bond_gam);
}

uni10::UniTensor<double> MPSInf::initLambda(int chi) {
  /// initialize a gamma tensor
  std::vector<uni10::Bond> bond_lam;
  bond_lam.push_back( uni10::Bond(uni10::BD_IN, chi) );
  bond_lam.push_back( uni10::Bond(uni10::BD_OUT, chi) );
  return uni10::UniTensor<double>(bond_lam);
}

void MPSInf::init() {
  ///
  if (gamma.size() > 0)
    gamma.clear();
  if (lambda.size() > 0)
    lambda.clear();

  int chi1, chi2;
  for (int i = 0; i < lat_size; ++i) {
    if (i%2 == 0) {
      chi1 = 1;
      chi2 = dim_phys;
    }
    else {
      chi1 = dim_phys;
      chi2 = 1;
    }
    gamma.push_back( initGamma(chi1, chi2, dim_phys) );
    lambda.push_back( initLambda(chi1) );
  }
}

void MPSInf::randomize() {
  /// randomize a complex MPS having only real part
  MPSInf::init();
  std::srand( time(NULL) );

  for (int i = 0; i < lat_size; ++i) {
    gamma[i].randomize();
    lambda[i].identity();
    lambda[i] *= (1./norm(lambda[i].getBlock()));
  }
}

uni10::UniTensor<double> MPSInf::expValAvg(uni10::UniTensor<double> op) {
  /// return expectation value (a contracted uni10) of an 1-site/2-site operator
  std::vector<uni10::UniTensor<double> > expV;
  expV.push_back( expVal(netLGLGL(lambda[0], gamma[0], lambda[1], gamma[1], lambda[0]), op) );
  expV.push_back( expVal(netLGLGL(lambda[1], gamma[1], lambda[0], gamma[0], lambda[1]), op) );
  return 0.5 * (expV[0] + expV[1]);
}

void MPSInf::mps2SiteSVD( uni10::UniTensor<double>& theta,
    uni10::UniTensor<double>& lam0, uni10::UniTensor<double>& gam0,
    uni10::UniTensor<double>& lam1, uni10::UniTensor<double>& gam1, uni10::UniTensor<double>& lam2,
    bool show_err ) {
  ///
  theta = permute( theta, theta.label(), theta.bondNum()/2 );
  std::vector<uni10::Matrix<double> > usv = svd( theta.getBlock() );
  double svsq_sum = norm(usv[1]);

  int dim_l = gam0.bond()[0].dim();
  int dim_m = std::min( (int)usv[1].col(), chi_max );
  int dim_r = gam1.bond()[2].dim();

  if ( dim_m != gam0.bond()[2].dim() )
    gam0 = initGamma(dim_l, dim_m, dim_phys );
  if ( dim_m != lam1.bond()[0].dim() )
    lam1 = initLambda(dim_m);
  if ( dim_m != gam1.bond()[0].dim() )
    gam1 = initGamma(dim_m, dim_r, dim_phys );

  uni10::UniTensor<double> l0i = tenInv(lam0);
  uni10::UniTensor<double> l2i = tenInv(lam2);

  resize( usv[0], dim_l * dim_phys, dim_m, uni10::INPLACE );
  resize( usv[1], dim_m, dim_m, uni10::INPLACE );
  resize( usv[2], dim_m, dim_phys * dim_r, uni10::INPLACE );

  gam0.putBlock( usv[0] );
  gam0 = netLG( l0i, gam0 );

  if (show_err) {
    double trunc_err = (svsq_sum - norm(usv[1])) / svsq_sum;
    std::cout << std::setprecision(12) << trunc_err << '\t';
  }
  usv[1] *= ( 1.0 / norm(usv[1]) );
  lam1.putBlock( usv[1] );

  gam1 = permute( gam1, gam1.label(), 1 );
  gam1.putBlock( usv[2] );
  gam1 = netGL( gam1, l2i );
}

//void MPSInf::idmrg(std::vector<uni10::UniTensor<double> >& mpo, int max_N, int lanczos_max_iter, double tolerance) {
void MPSInf::idmrg(std::vector<uni10::UniTensor<double> >& mpo, const idmrg_paras& paras){
  ///
  uni10::UniTensor<double> psi;
  double eng = 0.0;

  mpo.insert( mpo.end()-1, mpo[1] );
  uni10::UniTensor<double> ham2s = contrMPOLR( mpo );	// for energy output
  std::cout << ham2s << '\n';	// print 2-ste Hamiltonian

  fprintf(stdout, "Updating the matrix product states: \n\n");
  progressbar(0, 0, paras.max_N, true);

  uni10_uint64 iter = 0;
  int iter_limit = paras.lanczos_max_iter;

  for (int st = 0; st < paras.max_N; ++st) {

    if (st > 0) {
      std::swap( gamma[0], gamma[1] );
      std::swap( lambda[0], lambda[1] );
    }

    psi = netLGLGL( lambda[0], gamma[0], lambda[1], gamma[1], lambda[0] );
    psi = permute( psi, psi.label(), psi.bondNum() );
    double Norm = norm( psi.getBlock() );
    psi *= (1./Norm);
    iter_limit = std::min(paras.lanczos_max_iter, (int)psi.elemNum()+1);

    if (st == 0)
      iter = LanczosEigh( ham2s, psi, eng, iter_limit, paras.tolerance );

    else
      iter = LanczosEigh( mpo, psi, eng, iter_limit, paras.tolerance );
    MPSInf::mps2SiteSVD( psi, lambda[0], gamma[0], lambda[1], gamma[1], lambda[0] );

    // update mpo_l mpo_r
    renormMPOL( mpo, netLG( lambda[0], gamma[0] ), (st == 0) );
    renormMPOR( mpo, netGL( gamma[1], lambda[0] ), (st == 0) );

    if(st%paras.measure_per_n_iter==0){
      eng = MPSInf::expValAvg( ham2s )[0];
      progressbar(st+paras.measure_per_n_iter, 0, paras.max_N);
      cout.precision(8);
      cout.setf(ios::fixed, ios::floatfield);
      cout <<  ", ge: " << eng  << "\r";
      std::cout.flush();
    }
  }

  fprintf(stdout, "\n\n");

}

