#include "Utils.hpp"
#include "boost/math/special_functions.hpp"
#include "mkl.h"
#include <vector>
#include <cmath>
using namespace std;
using namespace uni10;
using namespace boost::math;

double Utils::GetMax(UniTensor<double> &T){
    return abs(T.getElem()[cblas_idamax(T.elemNum(), T.getElem() , 1)]);
}


uni10::UniTensor<double> Utils::MakeLocal(const double &J1, const double &J2, const double &Beta){

    // Geometry :
    //             u
    //            _|_
    //        l--|   |--r
    //           |   |
    //            ---
    //             |
    //             d
    //  Tensor : left-BD : in  /  right-BD : out
    //
    //            --------
    //         l--|      |--r
    //         u--|      |--d
    //            --------


    vector<Bond> Bds(2,Bond(BD_IN,2));
    Bds.push_back(Bond(BD_OUT,2));
    Bds.push_back(Bond(BD_OUT,2));

    UniTensor<double> O(Bds);

    vector<double> Raw(16);
    double sigD1, sigD2;
    // i-j-k-l = L-U-R-D
    for(int l=0;l<2;l++)
        for(int u=0;u<2;u++)
            for(int r=0;r<2;r++)
                for(int d=0;d<2;d++){
                    sigD1 = (2.*l - 1)*(2.*r - 1);
                    sigD2 = sigD1*(2.*l - 1)*(2.*u - 1);
                    Raw[l*8+u*4+r*2+d] = (1.+(2.*l-1.)*(2.*r-1.)*(2.*u-1.)*(2.*d-1.))/2
                                        * exp(-J1*Beta*0.5*(2.*(l+u+r+d)-4.) - J2*Beta*(sigD1+sigD2));

                }
    //O.SetElem(Raw);
    O.setElem(Raw);
    return O;

}

UniTensor<double> Utils::Make_T(const double &Beta,const double &h,const unsigned int &Chi){

    unsigned int Chi2 = Chi*Chi;
    unsigned int Chi3 = Chi2*Chi;

    ///Set Raw elements
    Matrix<double> rawN(Chi2,Chi2);

    for(unsigned int l=0;l<Chi;l++)
        for(unsigned int u=0;u<Chi;u++)
            for(unsigned int r=0;r<Chi;r++)
                for(unsigned int d=0;d<Chi;d++)
                    rawN[l*Chi3 + u*Chi2 + r*Chi + d] = sqrt(cyl_bessel_i(l,Beta)*
                                                             cyl_bessel_i(u,Beta)*
                                                             cyl_bessel_i(r,Beta)*
                                                             cyl_bessel_i(d,Beta))*
                                                        cyl_bessel_i(l+u-r-d,Beta*h);



    //rawN *= 1./norm(rawN);


    ///Set Tensor T_lurd:
    vector<Bond> Bds(2,Bond(BD_IN,Chi));
    Bds.push_back(Bond(BD_OUT,Chi));
    Bds.push_back(Bond(BD_OUT,Chi));

    UniTensor<double> T(Bds);
    //T.PutBlock(rawN);
    T.putBlock(rawN);
    return T;


}



void Utils::truncateLUs(const int dir, const int &chi, vector<UniTensor<double> >& svdUs,UniTensor<double> &T2){


	vector<int> ori_labels;
	vector<Bond> new_bonds;// = Us[0].bond();

	Matrix<double> blk;


	//L or U
	if(dir==0){
	    ori_labels = svdUs[1].label();
	    new_bonds = svdUs[1].bond();
	    new_bonds[1] = Bond(BD_OUT, chi);
	    resize(blk, svdUs[1].getBlock(), svdUs[1].getBlock().row(), chi, INPLACE);
	    svdUs[1].assign(new_bonds);
	    svdUs[1].putBlock(blk);
	    svdUs[1].setLabel(ori_labels);
	    T2 = contract(T2,svdUs[1],INPLACE);
        //cout << T2 ;
	    ori_labels = svdUs[3].label();
        svdUs[1].setLabel(ori_labels);
        T2 = contract(T2,svdUs[1],INPLACE);
        //cout << T2;
    }else{
        ori_labels = svdUs[3].label();
	    new_bonds = svdUs[3].bond();
	    new_bonds[1] = Bond(BD_OUT, chi);
	    resize(blk, svdUs[3].getBlock(), svdUs[3].getBlock().row(), chi, INPLACE);
	    svdUs[3].assign(new_bonds);
	    svdUs[3].putBlock(blk);
	    svdUs[3].setLabel(ori_labels);
	    T2 = contract(T2,svdUs[3],INPLACE);
        //cout << T2;
	    ori_labels = svdUs[1].label();
        svdUs[3].setLabel(ori_labels);
        T2 = contract(T2,svdUs[3],INPLACE);
        //cout << T2;
    }

}


//void Utils::Update(const int dir,const unsigned int &chi,UniTensor<double> &T, Network &Nwrk){
void Utils::Update(const int dir,const unsigned int &chi,UniTensor<double> &T, Network<double> &Nwrk){

	vector<int> per_lbl;
	if(dir == 0) per_lbl = vector<int>{0,1,3,4};
	else         per_lbl = vector<int>{4,0,1,3};

	vector<int> groups(4,1);

	vector<Matrix<double> > svdLs;
	vector<UniTensor<double> > svdUs;
	UniTensor<double> Core,T2;
    double e1=0,e2=0;
	//ContractArgs(T2,Nwrk,T,T);
	//T2.CombineBond({1,2});
	//T2.CombineBond({4,5});

	contract_args(T2,Nwrk,T,T);
    //cout << T2 << endl;


	T2.combineBond({1,2});
	T2.combineBond({4,5});


	//Hosvd( T2, T2.label(), groups, svdUs, Core, svdLs, INPLACE);

    //cout << "OK";
    //exit(1);

	///truncation :
	if(T2.bond(1).dim() > chi){
	    hosvd( T2, T2.label(), groups, svdUs, Core, svdLs, INPLACE);
        //cout << "IN" << endl;
        //cout << svdLs[1] << endl;
        //exit(1);
        for(unsigned int s=chi;s<svdUs[1].bond(1).dim();s++)
            e1 += pow(svdLs[1][s],2);
        for(unsigned int s=chi;s<svdUs[3].bond(1).dim();s++)
            e2 += pow(svdLs[3][s],2);
        //cout << e1 << " " << e2;
        //exit(1);
        if(e1<e2)
		    truncateLUs(0,chi,svdUs,T2);
        else
            truncateLUs(1,chi,svdUs,T2);

        T2.setLabel({0,3,1,4});
	}
    
    
	permute(T2,per_lbl,2,INPLACE);
    T.assign(T2.bond());
	T.putBlock(T2.getBlock());


}

