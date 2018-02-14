#ifndef __UNI10_ENV_SCALAPCK_MPI_H__
#define __UNI10_ENV_SCALAPCK_MPI_H__ 

#include <stdio.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_sys_info/uni10_sys_info_scalapack_mpi.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_blacs_wrapper_mpi.h"

namespace uni10{

  class uni10_env_mpi{

    public:

      uni10_env_mpi(): etype(mpi){};
      
      ~uni10_env_mpi();

      void clear();

      void set_communicate(uni10_exu_type _communicate);

      void used_memsize(const uni10_uint64& memsize);

      bool load_uni10_rc(int& _rank_mpi, int& _nprow_mpi, int& _npcol_mpi, int& _block_grid, int& _master);

      const sysinfo_mpi& get_info() const;

      friend void uni10_create(int argc, char** argv);

      friend void uni10_print_env_info();

    private:

      uni10_exu_type     etype;

      sysinfo_mpi        uni10_sys_info;

      void init(int argc=0, char** argv=NULL);

  };

  extern uni10_env_mpi env_variables;

};

#endif
