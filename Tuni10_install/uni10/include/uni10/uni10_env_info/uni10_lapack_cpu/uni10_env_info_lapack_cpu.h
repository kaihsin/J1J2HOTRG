#ifndef __UNI10_ENV_LAPACK_CPU_H__
#define __UNI10_ENV_LAPACK_CPU_H__ 

#include <stdio.h>

//#include "uni10/uni10_type.h"
//#include "uni10/uni10_error.h"
#include "uni10/uni10_sys_info/uni10_lapack_cpu/uni10_sys_info_lapack_cpu.h"

namespace uni10{

  class uni10_env_cpu{

    public:

      uni10_env_cpu(): etype(cpu){};

      ~uni10_env_cpu(){};

      void clear();

      void used_memsize(const uni10_uint64& memsize);

      bool load_uni10_rc();

      void print_env_info() const;

      const sysinfo_cpu& get_info() const;

      friend void uni10_create(int argc, char** argv);

      friend void uni10_print_env_info();
      
    private:

      uni10_exu_type    etype;

      sysinfo_cpu  uni10_sys_info;

      void init(int argc=0, char** argv=NULL);

  };

  extern uni10_env_cpu env_variables;

};

#endif
