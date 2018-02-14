#ifndef __UNI10_NETORDER_DEV_H__
#define __UNI10_NETORDER_DEV_H__

#include <utility>

#include "uni10/uni10_api/UniTensor.h"
#include "uni10/uni10_api/network_tools_dev/PseudoTensor.h"

typedef std::vector< std::pair<const void*, int> >  ary1d_uptr;
typedef std::vector< std::vector<uni10_int> > ary2d_label;
typedef std::vector< std::string >            ary1d_name;
typedef std::vector< uni10_int >              ary1d_order;

typedef std::vector< uni10::PseudoTensor_dev > ary1d_pten;
typedef std::vector< std::vector< uni10::PseudoTensor_dev> > ary2d_pten;

namespace uni10{

  class NetOrder_dev{

    public:

      NetOrder_dev();

      ~NetOrder_dev();

      NetOrder_dev(const ary1d_uptr& tens_type, const ary2d_label& label_arr, const ary1d_name& _names);

      ary1d_order generate_order();

      friend class PseudoTensor_dev;

    private:

      const ary1d_name* names;

      bool is_disjoint(const PseudoTensor_dev& T1, const PseudoTensor_dev& T2);

      bool is_overlap(const PseudoTensor_dev& T1, const PseudoTensor_dev& T2);

      PseudoTensor_dev psesudocontract(const PseudoTensor_dev& T1, const PseudoTensor_dev& T2);

      float get_cost(const PseudoTensor_dev& T1, const PseudoTensor_dev& T2);

      ary2d_pten tensor_set;

      float xi_min;

      int numT;

      std::vector<int> netorder_idx;

  };

};
#endif
