#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <cmath>
#include "uni10.hpp"
#include "Parser.hpp"
#include "Utils.hpp"
using namespace std;
using namespace uni10;

int main(int argc, char* argv[]){

    ///check arg :
    if(argc < 2){
        cout << "exec <.rc> \n";
        exit(1);
    }


    ///prepare parameters & parser
    unsigned int nL;
    double J1,J2,Beta;
    unsigned int chi;
    string ID;

    Parser pars;
    pars.Bind("ID",ID);
    pars.Bind("nL",nL);
    pars.Bind("J1",J1);
    pars.Bind("J2",J2);
    pars.Bind("Beta",Beta);
    pars.Bind("chi",chi);

    pars.Parse(argv[1]);
    pars.Check_All();
    pars.PrintVars();


    ///Prepare Network;
	Network Nwrk_lr("lr.net");
    Network Nwrk_ud("ud.net");

    ///Local Tensor:
    UniTensor<double> T = Utils::MakeLocal(J1,J2,Beta);


	///Cgran.
	//for(unsigned int itr=0;itr<nL;itr++){
		//printf("RC iter : %4d | L : %4d ",itr,(unsigned int)pow(2,itr+1) );
		Utils::Update(0,chi,T,Nwrk_lr);
		//Utils::Update(1,chi,T,Nwrk_ud);
	//}








    return 0 ;


}
