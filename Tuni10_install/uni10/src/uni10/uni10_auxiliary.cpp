#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_auxiliary.h"

namespace uni10{

  //env_type(env_info, _type) env_variables;

  void uni10_create(int argc, char** argv){

    env_variables.init(argc, argv);

  }

  void uni10_destroy(){

    env_variables.clear();

  }

  void uni10_print_env_info(){

    env_variables.print_env_info();

  }


}; // End of namespace
