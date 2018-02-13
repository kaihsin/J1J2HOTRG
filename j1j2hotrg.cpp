#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
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


    ///Prepare Buffer tensors;
    vector<UniTensor<double>> Us;
    UniTensor<double> S;
    vector<Matrix<double> > Ls;


    ///Local Tensor:
    UniTensor<double> T = Utils::MakeLocal(J1,J2,Beta);

    //Coursgran :

    //contract
    Network net("cg_lr.net");
    net.PutTensor(0,T);
    net.PutTensor(1,T);
    net.Launch(T);
    //hosvd
    Hosvd(T,vector<int>{1,4,3,6,0,5},vector<int>{2,2},Us,S,Ls, INPLACE);
    //truncate:









    return 0 ;


}
