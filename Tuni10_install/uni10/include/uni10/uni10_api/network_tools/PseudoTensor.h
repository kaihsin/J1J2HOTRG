#ifndef __UNI10_PSEUDOTENSOR_H__
#define __UNI10_PSEUDOTENSOR_H__

#include "uni10/uni10_api/UniTensor.h"

namespace uni10{

  class NetOrder;

  class PseudoTensor{

    public:
      PseudoTensor();

      PseudoTensor(std::vector<std::string>& _order, std::map<uni10_int32, uni10_int32>& _label_dim, 
          uni10_uint64 _bit = 0, uni10_float32 _cost = 0.0, uni10_bool _is_new = true);

      void printPseudoTensor() const;

      void netorder(char** netorder);

      void bufferInfo(uni10_int32* num_of_sub_tree, uni10_int32** size_of_sub_tree, char*** Tnames, uni10_int32 *numT);

      friend class NetOrder;

    private:

      std::vector<std::string> order;

      uni10_uint64 bit; 

      uni10_float32 cost;

      uni10_bool is_new;

      std::map<uni10_int32, uni10_int32> label_dim;

      uni10_int32 max_label;

  };

};

#endif
