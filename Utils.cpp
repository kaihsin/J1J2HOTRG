#include "Utils.hpp"
#include <vector>
#include <cmath>
using namespace std;
using namespace uni10;
uni10::UniTensor<double> Utils::MakeLocal(const double &J1, const double &J2, const double &Beta){

    // Geometry :
    //             j
    //            _|_
    //        i--|   |--k
    //           |   |
    //            ---
    //             |
    //             l
    //  Tensor : left-BD : in  /  right-BD : out
    //
    //            --------
    //         i--|      |--k
    //         j--|      |--l
    //            --------


    vector<Bond> Bds(2,Bond(BD_IN,2));
    Bds.push_back(Bond(BD_OUT,2));
    Bds.push_back(Bond(BD_OUT,2));

    UniTensor<double> O(Bds);

    vector<double> Raw(16);
    double sigD1, sigD2;
    // i-j-k-l = L-U-R-D
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            for(int k=0;k<2;k++)
                for(int l=0;l<2;l++){
                    sigD1 = (2.*i - 1)*(2.*k - 1);
                    sigD2 = sigD1*(2.*i - 1)*(2.*j - 1);
                    Raw[i*8+j*4+k*2+l] = (1.+(2.*i-1.)*(2.*k-1.)*(2.*j-1.)*(2.*l-1.))/2
                                        * exp(J1*Beta*0.5*(2.*(i+j+k+j)-4.) + J2*Beta*(sigD1+sigD2));

                }
    O.SetElem(Raw);
    return O;






}
