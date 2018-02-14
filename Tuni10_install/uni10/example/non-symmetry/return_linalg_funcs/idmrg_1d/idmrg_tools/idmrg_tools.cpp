#include "idmrg_tools.h"

uni10_uint64 LanczosEigh(
    std::vector<uni10::UniTensor<double> >& mpoH,
    uni10::UniTensor<double>& psi,
    double& E0, int max_iter, double err_tol );

uni10_uint64 LanczosEigh(
    uni10::UniTensor<double>& op,
    uni10::UniTensor<double>& psi,
    double& E0, int max_iter, double err_tol );

/*========================================

  Network contractions and block operations

  ========================================*/

uni10::UniTensor<double> netLG( uni10::UniTensor<double> la, uni10::UniTensor<double> ga ) {
  ///
  int lab_la[] = {0, 1};
  int lab_ga[] = {1, 100, 2};
  la.setLabel( lab_la );
  ga.setLabel( lab_ga );
  uni10::UniTensor<double> net = contract( la, ga, false );
  net = permute( net, net.label(), 2 );
  return net;
}

uni10::UniTensor<double> netGL( uni10::UniTensor<double> ga, uni10::UniTensor<double> la ) {
  ///
  int lab_ga[] = {0, 100, 1};
  int lab_la[] = {1, 2};
  ga.setLabel( lab_ga );
  la.setLabel( lab_la );
  uni10::UniTensor<double> net = contract( ga, la, false );
  return net;
}

uni10::UniTensor<double> netLGLGL(
    uni10::UniTensor<double> ll, uni10::UniTensor<double> ga,
    uni10::UniTensor<double> la, uni10::UniTensor<double> gb, uni10::UniTensor<double> lb ) {
  ///
  int lab_ll[] = {0, 1};
  int lab_ga[] = {1, 100, 2};
  int lab_la[] = {2, 3};
  int lab_gb[] = {3, 101, 4};
  int lab_lb[] = {4, 5};

  ll.setLabel( lab_ll );
  ga.setLabel( lab_ga );
  la.setLabel( lab_la );
  gb.setLabel( lab_gb );
  lb.setLabel( lab_lb );

  uni10::UniTensor<double> net = contract( ll, ga, false );
  net = contract( net, la, false );
  net = contract( net, gb, false );
  net = contract( net, lb, false );
  net = permute( net, net.label(), 3 );
  return net;
}

uni10::UniTensor<double> theta( uni10::UniTensor<double> ket, uni10::UniTensor<double> op ) {
  ///
  int lab_3b[] = {0, 100, 2};
  int lab_4b[] = {0, 100, 101, 5};
  int lab_op1[] = {200, 100};
  int lab_op2[] = {200, 201, 100, 101};
  int lab_net_3b[] = {0, 200, 2};
  int lab_net_4b1[] = {0, 200, 101, 5};
  int lab_net_4b2[] = {0, 200, 201, 5};

  uni10::UniTensor<double> net;

  if (op.bondNum() == 4 && ket.bondNum() == 4) {
    op.setLabel( lab_op2 );
    ket.setLabel( lab_4b );
    net = contract( op, ket, false );
    net = permute( net, lab_net_4b2, 3 );
  }
  else if (op.bondNum() == 2) {
    if (ket.bondNum() == 3) {
      op.setLabel( lab_op1 );
      ket.setLabel( lab_3b );
      net = contract( op, ket, false );
      net = permute( net, lab_net_3b, 2 );
    }
    else if (ket.bondNum() == 4) {
      op.setLabel( lab_op1 );
      ket.setLabel( lab_4b );
      net = contract( op, ket, false );
      net = permute( net, lab_net_4b1, 3 );
    }
  }
  return net;
}

uni10::UniTensor<double> expVal( uni10::UniTensor<double> ket, uni10::UniTensor<double> op ) {
  ///
  uni10::UniTensor<double> bra = ket;
  uni10::UniTensor<double> net = theta( ket, op );
  bra.setLabel( net.label() );
  bra = dagger( bra );
  net = contract( bra, net, false );
  return net;
}

uni10::UniTensor<double> tenInv( uni10::UniTensor<double>& ten ) {
  ///
  uni10::UniTensor<double> inv_ten = uni10::UniTensor<double>( ten.bond() );
  uni10::Matrix<double> blk = ten.getBlock();
  inv_ten.putBlock( inverse(blk) );
  return inv_ten;
}

uni10::UniTensor<double> mpoMatVec(
    std::vector<uni10::UniTensor<double> >& mpoH,
    uni10::UniTensor<double>& psi ) {
  ///
  uni10::UniTensor<double> net;

  int size = 4;
  int lab_l[] = {100, 0, 101};
  int lab_r[] = {size-2, size*100, size*100+1};

  std::vector<int> lab_m(4);
  std::vector<int> lab_psi, lab_net;
  for (int i = 1; i <= size; ++i) {
    lab_psi.push_back(i*100+1);
    lab_net.push_back(i*100);
  }
  psi.setLabel( lab_psi );
  mpoH[0].setLabel( lab_l );
  net = contract( mpoH[0], psi, false );

  for (int i = 2; i < size; ++i) {
    lab_m[0] = i-2;
    lab_m[1] = i*100;
    lab_m[2] = i-1;
    lab_m[3] = i*100+1;
    mpoH[i-1].setLabel( lab_m );
    net = contract( net, mpoH[i-1], false );
  }

  mpoH[size-1].setLabel( lab_r );
  net = contract( net, mpoH[size-1], false );
  net = permute(net, lab_net, size);

  return net;
}

uni10::UniTensor<double> contrMPOLR(
    std::vector<uni10::UniTensor<double> >& mpo ) {
  ///
  int lab_l[] = {100, 0, 101};
  int lab_r[] = {0, 200, 201};
  int lab_fin[] = {100, 200, 101, 201};

  mpo[0].setLabel( lab_l );
  mpo[mpo.size()-1].setLabel( lab_r );
  uni10::UniTensor<double> op = contract(mpo[0], mpo[mpo.size()-1], false);
  op = permute(op, lab_fin, 2);
  return op;
}

void renormMPOL(
    std::vector<uni10::UniTensor<double> >& mpoH,
    uni10::UniTensor<double> ket, bool initial ) {
  ///
  uni10::UniTensor<double> bra = dagger(ket);

  int lab_new[] = {0, 1, 2};
  int lab_npm[] = {0, 1, 2, -1};

  if (initial) {
    int lab_bra[] = {0, -10, 200};
    int lab_l[]   = {200, 1, 100};
    int lab_ket[] = {-10, 100, 2};
    bra.setLabel( lab_bra );
    ket.setLabel( lab_ket );

    mpoH[0].setLabel( lab_l );
    mpoH[0] = contract( bra, mpoH[0], false );
    mpoH[0] = contract( mpoH[0], ket, false );
    mpoH[0] = permute( mpoH[0], lab_new, 1 );
  }
  else {
    int lab_bra[] = {0, 200, 201};
    int lab_l[]   = {200, 10, 100};
    int lab_m[]   = {10, 201, 1, 101};
    int lab_ket[] = {100, 101, 2};
    bra.setLabel( lab_bra );
    ket.setLabel( lab_ket );

    uni10::UniTensor<double> lvec;
    lvec = mpoH[0];
    lvec.setLabel( lab_l );
    mpoH[1].setLabel( lab_m );

    mpoH[0] = contract( bra, lvec, false );
    mpoH[0] = contract( mpoH[0], mpoH[1], false );
    mpoH[0] = contract( mpoH[0], ket, false );
    mpoH[0] = permute( mpoH[0], lab_new, 1 );
  }
}

void renormMPOR(
    std::vector<uni10::UniTensor<double> >& mpoH,
    uni10::UniTensor<double> ket, bool initial ) {
  ///
  uni10::UniTensor<double> bra = dagger(ket);

  int n = (int)mpoH.size()-1;
  int lab_new[] = {1, 0, 2};

  if (initial) {
    int lab_bra[] = {-10, 0, 200};
    int lab_r[]   = {1, 200, 100};
    int lab_ket[] = {2, 100, -10};
    bra.setLabel( lab_bra );
    ket.setLabel( lab_ket );

    mpoH[n].setLabel( lab_r );
    mpoH[n] = contract( bra, mpoH[n], false );
    mpoH[n] = contract( mpoH[n], ket, false );
    mpoH[n] = permute( mpoH[n], lab_new, 2 );
  }
  else {
    int lab_bra[] = {200, 0, 201};
    int lab_r[]   = {10, 200, 100};
    int lab_m[]   = {1, 201, 10, 101};
    int lab_ket[] = {2, 101, 100};
    bra.setLabel( lab_bra );
    ket.setLabel( lab_ket );

    uni10::UniTensor<double> rvec;
    rvec = mpoH[n];
    rvec.setLabel( lab_r );
    mpoH[n-1].setLabel( lab_m );

    mpoH[n] = contract( bra, rvec, false );
    mpoH[n] = contract( mpoH[n], mpoH[n-1], false );
    mpoH[n] = contract( mpoH[n], ket, false );
    mpoH[n] = permute( mpoH[n], lab_new, 2 );
  }
}

/*========================================

  Lanczos method

  ========================================*/

uni10_uint64 LanczosEigh(
    std::vector<uni10::UniTensor<double> >& mpoH, uni10::UniTensor<double>& psi,
    double& E0, int max_iter, double err_tol ) {
  /// psi: all in-bonds
  uni10::Matrix<double> psi_blk = psi.getBlock();

  uni10_lanczos_custom3_paras<double> cust3_paras(mpoH, psi, psi_blk);
  uni10_lanczos_paras<double, uni10_lanczos_custom3_paras<double> > l_paras(cust3_paras);

  int info;
  uni10_uint64 n=0;
  E0 = Lanczos(l_paras, n, max_iter, 2, info, err_tol);
  psi.putBlock(psi_blk);
  return n;
}

uni10_uint64 LanczosEigh(
    uni10::UniTensor<double>& op, uni10::UniTensor<double>& psi,
    double& E0, int max_iter, double err_tol ) {
  /// psi: all in-bonds
  uni10::Matrix<double> op_blk = op.getBlock();
  uni10::Matrix<double> psi_blk = psi.getBlock();

  uni10_lanczos_default_paras<double> default_paras(op_blk, psi_blk);
  uni10_lanczos_paras<double, uni10_lanczos_default_paras<double> > l_paras(default_paras);

  int info;
  uni10_uint64 n=0;
  E0 = Lanczos(l_paras, n, max_iter, 2, info, err_tol);
  psi.putBlock(psi_blk);
  return n;
}

//======================================

