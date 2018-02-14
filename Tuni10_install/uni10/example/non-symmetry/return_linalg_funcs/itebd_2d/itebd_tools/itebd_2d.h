#ifndef __ITEBD_2D_H__
#define __ITEBD_2D_H__

#include "uni10/uni10.hpp"
#include "../../../../common/common_tools.h"

using namespace std;
using namespace uni10;

struct itebd_paras{

  itebd_paras(){
    D = 5; max_D = 0;
    max_N = 1000000;
    tau = 1.0e-4; eps = 1.0e-14; cut_off = -1.;
  }

  void load_itebd_paras(){

    FILE* rcfp = fopen(".itebdrc", "r");

    int max_len = 256;
    char buffer[max_len];


    char* pch;
    while(fgets(buffer, max_len, rcfp)){

      pch = strtok(buffer, ":");
      pch = strtok(NULL, ":");

      if(strcmp ("D", buffer)==0)
        D = atoi(pch);

      else if(strcmp ("max_D", buffer)==0)
        max_D = atoi(pch);

      else if(strcmp ("max_N", buffer)==0)
        max_N = atoi(pch);

      else if(strcmp ("tau", buffer)==0)
        tau = atof(pch);

      else if(strcmp ("eps", buffer)==0)
        eps = atof(pch);

      else if(strcmp ("cut_off", buffer)==0)
        cut_off = atof(pch);

      else if(strcmp ("measure_per_n_iter", buffer)==0)
        measure_per_n_iter = atoi(pch);

      else if(buffer[0] =='#' || pch == NULL)
        continue;

      else{
        fprintf(stdout, "%s", "Setting the parameters with wrong names.");
        exit(0);
      }

    }

    fclose(rcfp);

  }

  void print_info() const{

    fprintf(stdout, "=====================================\n");
    fprintf(stdout, "|        Parameters of iTEBD        |\n");
    fprintf(stdout, "=====================================\n");
    fprintf(stdout, "# D      : %d\n", D);
    fprintf(stdout, "# max_D  : %d\n", max_D);
    fprintf(stdout, "# max_N  : %lld\n", max_N);
    fprintf(stdout, "# tau    : %.5f\n", tau);

    if(eps > 0)
      fprintf(stdout, "# eps    : 1e%.1f\n", log10(eps));
    else
      fprintf(stdout, "# eps    : %.2f\n", eps);

    if(cut_off > 0)
      fprintf(stdout, "# cut_off: 1e%.1f\n", log10(cut_off));
    else
      fprintf(stdout, "# cut_off: %.2f\n", cut_off);

    fprintf(stdout, "# measure_per_n_iter: %d\n", measure_per_n_iter);
    fprintf(stdout, "=====================================\n");
  }

  int D, max_D, measure_per_n_iter;     // D: Virtual bonds dimension.
                       // d: Physical bonds dimension.
                       // max_D: Upper bondary of vitual bonds dimesnion.
  long long int max_N; // max_N: Upper bondary of update interations.
  double tau, eps, cut_off;

};

template<typename T>
class iTEBD_2D{

  public:
    iTEBD_2D(const UniTensor<T>& _H, const itebd_paras& paras);

    ~iTEBD_2D();

    void setHamiltonian(const UniTensor<T>& _H);

    UniTensor<T> get_gate(const UniTensor<T>& H);

    void Optimize();
        

  private:
    uni10_int dim;      // Physical bond dimension.
    uni10_int D;        // Virtual bond dimension.
    uni10_double64 tau; // trotter constant
    uni10_double64 eps; // 
    uni10_uint64 max_N; // Maximun iteration number.
    uni10_int cut_off;  // The cut off for dynamic truncation.
    uni10_int max_D;    // Maximum virtual bond dimension after truncation.
    uni10_int measure_per_n_iter;

    vector<UniTensor<T> > gammas;
    vector<Matrix<T> > lambdas;
	
    UniTensor<T> H;       // Hamiltonian filled of real components.

    void init();

    void bondcatAll( const int idir );
    void bondrmAll( const int idir );
    UniTensor<T> contractGammas( const int idir );
    void contractTwoSiteOp( UniTensor<T> &theta, UniTensor<T> &twoSiteOp );
    void updateGammas( Matrix<T> &u, Matrix<T> &vT, const int idir );
    T measureNorm( const int idir );
    T measureExpe( UniTensor<T> &twoSiteOp, const int idir );
};

#endif
