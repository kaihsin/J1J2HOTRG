#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <cmath>
#include <fstream>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
//#include "uni10.hpp"
#include "uni10/uni10.hpp"
#include "Parser.hpp"
#include "Utils.hpp"
using namespace std;
using namespace uni10;

#define CheckSeg(cond,errmsg) if(cond){printf("%s\n",errmsg); exit(1);}
struct Param{
    double J1,J2,T;
};

int main(int argc, char* argv[]){

    ///check arg :
    if(argc < 3){
        cout << "exec <.rc> <Tls>\n";
        exit(1);
    }


    ///prepare parameters & parser
    unsigned int nL;
    unsigned int chi;
    string ID;
    vector<Param> parls;

    Parser pars;
    pars.Bind("ID",ID);
    pars.Bind("nL",nL);
    pars.Bind("chi",chi);

    pars.Parse(argv[1]);
    pars.Check_All();
    pars.PrintVars();

    ///Read Parameters:
    fstream f;
    double tmpT,tmpJ1,tmpJ2;
    f.open(argv[2],ios::in);
    CheckSeg(!f.is_open(),"[ERROR] open Tls fail.");
    while(f>>tmpJ1>>tmpJ2>>tmpT){
        Param pars = {tmpJ1,tmpJ2,tmpT};
        parls.push_back(pars);
    }

    ///Prepare Network;
	//Network Nwrk_lr("lr.net");
    //Network Nwrk_ud("ud.net");
    Network<double> Nwrk_lr("lr.net");
    Network<double> Nwrk_ud("ud.net");

    ///SaveDir:
    string datDir = "Data/" + ID;
    mkdir("Data",S_IRWXU);
    mkdir(datDir.c_str(),S_IRWXU);
    string savpath = datDir + "/data";
    ofstream cfout;
    cfout.open(savpath.c_str(),ios::out| ios::trunc);
    cfout.close();

    ///buffers:
    double nrm;
    vector<double> lnNrms;
    unsigned int N;
    double tmpZ;
    double lnZ;

    for(unsigned int p=0;p<parls.size();p++){

        printf("J1 = %05.8lf \t J2 = %05.8lf \t T = %05.8lf \n",parls[p].J1,parls[p].J2,parls[p].T);

        ///Local Tensor :
        UniTensor<double> T = Utils::MakeLocal(parls[p].J1,parls[p].J2,(double)1./parls[p].T);


        //normalize:
        //nrm = Norm(T.GetBlock());
        //nrm = trace(T.getBlock());
        nrm = abs(Utils::GetMax(T));
        lnNrms.push_back(log(nrm));
        T*= (double)1./nrm;

        ///Cgran.
        for(unsigned int itr=0;itr<nL;itr++){
            //printf(" cgram : %4d | L : %4d \n",itr,(unsigned int)pow(2,itr+1) );
            Utils::Update(0,chi,T,Nwrk_lr);
            Utils::Update(1,chi,T,Nwrk_ud);
            nrm = norm(T.getBlock());
            //nrm = Utils::GetMax(T);
            //nrm = abs(trace(T.getBlock()));
            cout << "nrm " << nrm << endl;
            //nrm /= 2;
			T *= (double)1./nrm;
            //if(nrm>100)
            //    T *= (double)1./nrm;
            //else 
            //    nrm = 1;
            lnNrms.push_back(log(nrm));
        }

        ///Obsv:
        double L = pow(2,nL);

        ///Calc PTfx (Z)
        T = partialTrace(T,0,2);
        tmpZ = trace(T.getBlock());
        cout << tmpZ << endl;
        lnZ = log(tmpZ)/L/L;
        cout << lnZ << endl;
        unsigned int itr = 0;
        while(lnNrms.size()){
            lnZ += lnNrms.back() * (pow(4,itr)/L/L);
            cout << " " << lnNrms.back() ;
            itr++;
            lnNrms.pop_back();
        }
        cout << endl;


        //lnZ = ln_nrmZ;
        //F   = -ln_nrmZ/Beta/(L*L);
        cfout.open(savpath.c_str(),ios::out| ios::app);
        cfout << fixed << setprecision(14) << parls[p].J1   << " "
                                           << parls[p].J2   << " "
                                           << parls[p].T    << " "
                                           << lnZ           << " "
                                           << -lnZ*parls[p].T << endl;
        cfout.close();
        //cfout << "%11.11lf %11.11lf %11.11lf"
        printf("lnZ = %010.8lf \t F(per site) = %3.8lf\n",lnZ,-lnZ*parls[p].T);

    }




    return 0 ;


}

