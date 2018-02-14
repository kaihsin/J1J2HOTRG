#ifndef MPSINF_HPP
#define MPSINF_HPP

#include <string>

#include "uni10/uni10.hpp"

#include "idmrg_1d.h"
#include "idmrg_tools.h"
#include "../../../../common/common_tools.h"

class MPSInf {

  public:
    /// constructor
    MPSInf(int X, int d, int L = 2);
    ~MPSInf();

    void init();
    void randomize();

    //void idmrg(std::vector<uni10::UniTensor<double>>& mpo, int max_N, int lanczos_max_iter, double tolerance);
    void idmrg(std::vector<uni10::UniTensor<double>>& mpo, const idmrg_paras& paras);

    uni10::UniTensor<double> expValAvg(uni10::UniTensor<double> op);

  private:
    int lat_size;
    int dim_phys;
    int chi_max;
    std::vector<uni10::UniTensor<double>> gamma;
    std::vector<uni10::UniTensor<double>> lambda;
    uni10::UniTensor<double> initGamma(int chi1, int chi2, int d);
    uni10::UniTensor<double> initLambda(int chi);

    void mps2SiteSVD( uni10::UniTensor<double>& theta,
        uni10::UniTensor<double>& lam0, uni10::UniTensor<double>& gam0,
        uni10::UniTensor<double>& lam1, uni10::UniTensor<double>& gam1, uni10::UniTensor<double>& lam2,
        bool show_err = false );

};


#endif
