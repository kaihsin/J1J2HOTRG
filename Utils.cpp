#include "Utils.hpp"
#include <vector>
#include <cmath>
using namespace std;
using namespace uni10;
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
                                        * exp(J1*Beta*0.5*(2.*(l+u+r+d)-4.) + J2*Beta*(sigD1+sigD2));

                }
    O.SetElem(Raw);
    return O;

}


void Utils::truncateLUs(const int dir, const int &chi, vector<UniTensor<double> >& svdUs,UniTensor<double> &T2){


	vector<int> ori_labels;
	vector<Bond> new_bonds;// = Us[0].bond();

	Matrix<double> blk;


	//L or U
	//cout << "svd1\n";
	//svdUs[1].printDiagram();
	ori_labels = svdUs[1].label();
	new_bonds = svdUs[1].bond();
	new_bonds[1] = Bond(BD_OUT, chi);
	Resize(blk, svdUs[1].GetBlock(), svdUs[1].GetBlock().row(), chi, INPLACE);
	svdUs[1].Assign(new_bonds);
	svdUs[1].PutBlock(blk);
	svdUs[1].SetLabel(ori_labels);
	Contract(T2,T2,svdUs[1],INPLACE);

	ori_labels = svdUs[3].label();
	new_bonds = svdUs[3].bond();
	new_bonds[1] = Bond(BD_OUT, chi);
	Resize(blk, svdUs[3].GetBlock(), svdUs[3].GetBlock().row(), chi, INPLACE);
	svdUs[3].Assign(new_bonds);
	svdUs[3].PutBlock(blk);
	svdUs[3].SetLabel(ori_labels);
	Contract(T2,T2,svdUs[3],INPLACE);


}


void Utils::Update(const int dir,const unsigned int &chi,UniTensor<double> &T, Network &Nwrk){

	vector<int> per_lbl;
	if(dir == 0) per_lbl = vector<int>{0,1,3,4};
	else         per_lbl = vector<int>{4,0,1,3};

	vector<int> groups(4,1);

	vector<Matrix<double> > svdLs;
	vector<UniTensor<double> > svdUs;
	UniTensor<double> Core,T2;

	ContractArgs(T2,Nwrk,T,T);
	T2.CombineBond({1,2});
	T2.CombineBond({4,5});
	//T2.printDiagram();
	Hosvd( T2, T2.label(), groups, svdUs, Core, svdLs, INPLACE);
	//exit(1);
	//cout << svdUs[1].bond(1).dim() << endl;
	///truncation :
	if(svdUs[1].bond(1).dim() > chi){
		truncateLUs(dir,chi,svdUs,T2);
	}else{
		Contract(T2,T2,svdUs[1],INPLACE);
		Contract(T2,T2,svdUs[3],INPLACE);
	}

	T2.SetLabel({0,3,1,4});

	Permute(T2,per_lbl,2,INPLACE);
	T.PutBlock(T2.GetBlock());
}

