#ifndef __IDMRG_1D_H__
#define __IDMRG_1D_H__

#include "uni10/uni10.hpp"
#include "../../../../common/common_tools.h"

using namespace std;
using namespace uni10;

struct idmrg_paras{

  idmrg_paras(){

    chi        = 5;
    max_N      = 100;
    tolerance  = 1.0e-15;
    lanczos_max_iter = 500;
    measure_per_n_iter = 1;

  }

  void load_idmrg_paras(){

    FILE* rcfp = fopen(".idmrgrc", "r");

    int max_len = 256;
    char buffer[max_len];

    char* pch;
    while(fgets(buffer, max_len, rcfp)){

      pch = strtok(buffer, ":");
      pch = strtok(NULL, ":");

      if(strcmp ("chi", buffer)==0)
        chi = atoi(pch);

      else if(strcmp ("max_N", buffer)==0)
        max_N = atoi(pch);

      else if(strcmp ("lanczos_max_iter", buffer)==0)
        lanczos_max_iter = atoi(pch);

      else if(strcmp ("tol", buffer)==0)
        tolerance = atof(pch);

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
    fprintf(stdout, "|        Parameters of iDMRG        |\n");
    fprintf(stdout, "=====================================\n");
    fprintf(stdout, "# chi        : %d\n"  , chi);
    fprintf(stdout, "# max_N      : %d\n"  , max_N);
    fprintf(stdout, "# tolerance  : %.5f\n", tolerance);
    fprintf(stdout, "# lanczos_max_iter   : %d\n", lanczos_max_iter);
    fprintf(stdout, "# measure_per_n_iter : %d\n", measure_per_n_iter);
    fprintf(stdout, "=====================================\n");
  }

  int chi;              // chi: Virtual bonds dimension.
  int max_N;            // max_N: Upper bound of DMRG steps.
  int lanczos_max_iter; // lanczos_max_iter: Upper bound of Lanczos interations in each DMRG step.
  double tolerance;
  int measure_per_n_iter;

};

#endif
