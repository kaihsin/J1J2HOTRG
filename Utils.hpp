#ifndef UTILS_HPP_INCLUDED
#define UTILS_HPP_INCLUDED

#include <iostream>
#include <cstdlib>
#include <cstring>
#include "uni10.hpp"


namespace Utils{

    uni10::UniTensor<double> MakeLocal(const double &J1, const double &J2, const double &Beta);
	void truncateLUs(const int dir, const int &chi, std::vector<uni10::UniTensor<double> >& svdUs,uni10::UniTensor<double> &T2);
    void Update(const int dir,const unsigned int &chi,uni10::UniTensor<double> &T, uni10::Network &Nwrk);



}




#endif // UTILS_HPP_INCLUDED
