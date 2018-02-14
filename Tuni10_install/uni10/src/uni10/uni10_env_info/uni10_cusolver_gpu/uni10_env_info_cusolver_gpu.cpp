#include "uni10/uni10_error.h"
#include "uni10/uni10_env_info/uni10_cusolver_gpu/uni10_cusolver_gpu_initDconst.h"
#include "uni10/uni10_env_info/uni10_cusolver_gpu/uni10_env_info_cusolver_gpu.h"

namespace uni10{

  uni10_env_gpu env_variables;

  void uni10_env_gpu::init(int argc, char** argv){
    
    uni10_sys_info.default_init();

    if(!this->load_uni10_rc(uni10_sys_info)){


    }

    uni10_init_device_const(this->uni10_sys_info);


  }

  void uni10_env_gpu::used_memsize(const uni10_uint64& memsize){

    uni10_sys_info.free_byte -= memsize;

  }

  void uni10_env_gpu::clear(){

    etype         = gpu;
    uni10_sys_info.clear();

  }

  bool uni10_env_gpu::default_ongpu_flag(){

    bool ongpu;

    if(uni10_sys_info.runtime_type == only_cpu){

      ongpu = false;

    }

    else if(uni10_sys_info.runtime_type == only_gpu){

      ongpu = true;

    }

    else if(uni10_sys_info.runtime_type == hybrid){

      uni10_error_msg(true, "%s", "Developing");
      ongpu = false;

    }

    return ongpu;

  }

  bool uni10_env_gpu::load_uni10_rc(sysinfo_gpu& _uni10_sys_info){

    bool exsist_rc = true;
    FILE* rcfp = NULL;

    // 1.Try to open the .uni10rc in home direction.
    fopen("~/.uni10rc", "r");
    // 2.Try to open the .uni10rc in the loacl direction.
    if(!rcfp)
      rcfp = fopen(".uni10rc", "r");

    if(!rcfp)
      exsist_rc = false;
    else{

      this->etype = gpu;

      int max_len = 256;
      char buffer[max_len];

      char* pch;
      while(fgets(buffer, max_len, rcfp)){

        pch = strtok(buffer, ":");
        pch = strtok(NULL, ":");

        if(strcmp ("MODE", buffer)==0){

          int mode = atoi(pch);
          if(mode == 0)
            _uni10_sys_info.runtime_type = only_cpu;

          else if(mode == 1)
            _uni10_sys_info.runtime_type = hybrid;

          else if(mode == 2)
            _uni10_sys_info.runtime_type = only_gpu;
        }

        else if(strcmp ("THREADSPERBLOCK_X", buffer)==0)
          _uni10_sys_info.threadsPerBlock_x = atoi(pch);

        else if(strcmp ("THREADSPERBLOCK_Y", buffer)==0)
          _uni10_sys_info.threadsPerBlock_y = atoi(pch);

        else
          uni10_error_msg(true, "%s", "Setting the parameters with wrong names.");

      }

    }

    return exsist_rc;

  }

  const sysinfo_gpu& uni10_env_gpu::get_info() const{

    return uni10_sys_info;

  }

  void uni10_env_gpu::print_env_info() const{

    uni10_sys_info.print_sys_info();
    std::map<std::string, uni10_int> dev_info = uni10_get_device_const();
    std::map<std::string, uni10_int>::const_iterator it = dev_info.begin();

    fprintf(stdout,"\n----- DEVICE CONSTANT -----\n");
    for(; it != dev_info.end(); it++)
      fprintf(stdout, "%s: %d\n", it->first.c_str(), it->second);
    fprintf(stdout,"---------------------------\n");

  }

}
