#include <limits.h>

#include "uni10/uni10_api/network_tools_dev/PseudoTensor.h"

namespace uni10{

  PseudoTensor_dev::PseudoTensor_dev(): bit(0), cost(0), is_new(true), max_label(INT_MIN){};

  PseudoTensor_dev::PseudoTensor_dev(std::vector<uni10_int>& _order_idx, std::map<uni10_int, uni10_int>& _label_dim, 
      uni10_uint64 _bit, uni10_float32 _cost, uni10_bool _is_new):
  order_idx(_order_idx), bit(_bit), cost(_cost), is_new(_is_new), label_dim(_label_dim){

    max_label = label_dim.rbegin()->first;

  }

  PseudoTensor_dev::~PseudoTensor_dev(){

  };

  void PseudoTensor_dev::printPseudoTensor(const std::vector<std::string>* correspond_name) const {

    fprintf(stdout, "\n");
    fprintf(stdout, "=======================\n");

    for(uni10_int i = 0; i < (uni10_int)order_idx.size(); i++){
      if(i == 0)
        fprintf(stdout, "ORDER_IDX: %d", order_idx[i]);
      else
        fprintf(stdout, " %d", order_idx[i]);
    }
    fprintf(stdout, "\n");

    if(correspond_name != NULL)
      for(uni10_int i = 0; i < (uni10_int)order_idx.size(); i++){
        if(i == 0)
          fprintf(stdout, "ORDER_IDX: %s", (*correspond_name)[order_idx[i]].c_str());
        else if(order_idx[i] == -1)
          fprintf(stdout, " #");
        else
          fprintf(stdout, " %s", (*correspond_name)[order_idx[i]].c_str());
      }

    fprintf(stdout, "\n");
    fprintf(stdout, "Bit : %ld\n", bit);
    fprintf(stdout, "Cost: %.7f\n", cost);
    if(is_new)
      fprintf(stdout, "New : True\n");
    else
      fprintf(stdout, "New : False\n");

    if(label_dim.size() != 0){

      fprintf(stdout, "{label, dim}: ");
      std::map<uni10_int, uni10_int>::const_iterator it = label_dim.begin();
      fprintf(stdout, "{%d, %d}", it->first, it->second);
      for(it++; it != label_dim.end(); it++)
        fprintf(stdout, ", {%d, %d}", it->first, it->second);
      fprintf(stdout, "\n");
      fprintf(stdout, "Max_label: %d\n", max_label);

    }else
      fprintf(stdout, "Scaler !!\n");
    fprintf(stdout, "=======================\n");
    fprintf(stdout, "\n");

  }

};
