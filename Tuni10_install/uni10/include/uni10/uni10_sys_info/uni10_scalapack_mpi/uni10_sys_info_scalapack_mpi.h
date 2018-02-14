#ifndef __UNI10_SYS_INFO_SCALAPACK_MPI_H__
#define __UNI10_SYS_INFO_SCALAPACK_MPI_H__

//#define OSX

#if defined(LINUX)
#include <sys/sysinfo.h>
#elif defined(OSX)
#include <sys/vmmeter.h>
#endif

#include <sys/types.h>
#include <sys/sysctl.h>

#include <unistd.h>

#include "mpi.h"

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"

namespace uni10{

  struct sysinfo_mpi{

    // Developping 
    sysinfo_mpi(): rank_mpi(0), nprocs_mpi(1), master(0), 
      ictxt(-1), nprow(1), npcol(1), blockgrid(32), myrow(0), mycol(0){};

    sysinfo_mpi(uni10_int& _rank_mpi, uni10_int& _nprocs_mpi, uni10_int& _master,
       uni10_int& _ictxt ,uni10_int& _nprow, uni10_int& _npcol, uni10_int& _blockgrid, uni10_int& _myrow, uni10_int& _mycol): 
      rank_mpi(_rank_mpi), nprocs_mpi(_nprocs_mpi), master(_master), 
      ictxt(_ictxt), nprow(_nprow), npcol(_npcol), blockgrid(_blockgrid), myrow(_myrow), mycol(_mycol){};

    void clear(){
      //rank_mpi    = 0; 
      //nprocs_mpi  = 2; 
      //blockgrid   = 16;
      //master      = 0;
      //ictxt       = -1;
      //nprow       = 1; 
      //npcol       = 1; 
      //myrow       = 0;
      //mycol       = 0;
    };

    void print_env_info(){

      if(rank_mpi==0){

        fprintf(stdout, "\n#######  Uni10 environment information  #######\n");
        fprintf(stdout, "# MPI   Processes : %d \n", nprocs_mpi);
        fprintf(stdout, "# Row   Processes : %d \n", nprow);
        fprintf(stdout, "# Col   Processes : %d \n", npcol);
        fprintf(stdout, "# Blocks Grid Size: %d \n", blockgrid);
        fprintf(stdout, "###############################################\n\n");

      }

      MPI_Barrier(MPI_COMM_WORLD);

    };

    //
    // Parameters
    //
    // For mpi.
    uni10_int rank_mpi, nprocs_mpi, master;

    // For blacs.
    uni10_int ictxt, nprow, npcol, blockgrid, myrow, mycol;
  };

}

#endif
