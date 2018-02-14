#ifndef __UNI10_ENV_CUSOLVER_GPU_H__
#define __UNI10_ENV_CUSOLVER_GPU_H__ 

#include <stdio.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_sys_info/uni10_cusolver_gpu/uni10_sys_info_cusolver_gpu.h"
#include "uni10/uni10_env_info/uni10_cusolver_gpu/uni10_cusolver_gpu_initDconst.h"

namespace uni10{

  class uni10_env_gpu{

    public:

      uni10_env_gpu(): etype(gpu){};

      ~uni10_env_gpu(){};

      void clear();

      void used_memsize(const uni10_uint64& memsize);

      bool default_ongpu_flag();

      bool load_uni10_rc( sysinfo_gpu& _uni10_sys_info );

      void print_env_info() const;

      const sysinfo_gpu& get_info() const;

      friend void uni10_create(int argc, char** argv);

      friend void uni10_print_env_info();

    private:

      uni10_exu_type             etype;

      sysinfo_gpu                uni10_sys_info;

      void init(int argc=0, char** argv=NULL);

  };

  extern uni10_env_gpu env_variables;

};

#endif
