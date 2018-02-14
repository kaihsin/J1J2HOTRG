#ifndef __UNI10_PSEUDOTENSOR_DEV_H__
#define __UNI10_PSEUDOTENSOR_DEV_H__

#include "uni10/uni10_api/UniTensor.h"

namespace uni10{

  class NetOrder_dev;

  class PseudoTensor_dev{

    public:

      PseudoTensor_dev();

      ~PseudoTensor_dev();

      PseudoTensor_dev(std::vector<uni10_int>& _order_idx, std::map<uni10_int, uni10_int>& _label_dim, 
          uni10_uint64 _bit = 0, uni10_float32 _cost = 0.0, uni10_bool _is_new = true);

      void printPseudoTensor(const std::vector<std::string>* correspond_names=NULL) const;

      friend class NetOrder_dev;

    private:

      std::vector<uni10_int> order_idx;

      uni10_uint64 bit; 

      uni10_float32 cost;

      uni10_bool is_new;

      std::map<uni10_int, uni10_int> label_dim;

      uni10_int max_label;

  };

};

#endif
