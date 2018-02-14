#include "uni10/uni10_error.h"
#include "uni10/uni10_env_info/uni10_lapack_cpu/uni10_env_info_lapack_cpu.h"

namespace uni10{

  uni10_env_cpu env_variables;

  void uni10_env_cpu::init(int argc, char** argv){
    
    argc = 0;
    argv = NULL;
    etype         = cpu;
    this->uni10_sys_info.init();
    this->load_uni10_rc();

  }

  void uni10_env_cpu::clear(){

    etype         = no;
    uni10_sys_info.clear();

  }

  void uni10_env_cpu::used_memsize(const uni10_uint64& memsize){

    this->uni10_sys_info.free_memsize -= memsize;

  }

  bool uni10_env_cpu::load_uni10_rc(){

    bool exsist_rc = true;
    FILE* rcfp = fopen("~/.uni10rc", "r");

    if(!rcfp)
      exsist_rc = false;
    else{

      this->etype = cpu;

      uni10_error_msg(true, "%s", "Developping !!!");

    }

    return exsist_rc;
    
  }

  const sysinfo_cpu& uni10_env_cpu::get_info() const{

    return uni10_sys_info;

  }

  void uni10_env_cpu::print_env_info() const{

    uni10_sys_info.print_sys_info();

  }

}
