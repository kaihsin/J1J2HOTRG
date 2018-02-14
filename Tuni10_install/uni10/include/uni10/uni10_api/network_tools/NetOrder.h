#ifndef __UNI10_NETORDER_H__
#define __UNI10_NETORDER_H__

#include "uni10/uni10_api/UniTensor.h"
#include "uni10/uni10_api/network_tools/PseudoTensor.h"

typedef std::vector< uni10::UniTensor< uni10_double64>* >   ary1d_uptr_d;
typedef std::vector< uni10::UniTensor< uni10_complex128>* > ary1d_uptr_z;
typedef std::vector< uni10::PseudoTensor > ary1d_pten;
typedef std::vector< std::vector< uni10::PseudoTensor> > ary2d_pten;

namespace uni10{

  class NetOrder{

    public:

      NetOrder();

      NetOrder(const ary1d_uptr_d& ts);

      NetOrder(const ary1d_uptr_z& ts);

      ~NetOrder();

      char* generate_order();

    private:

      bool is_disjoint(const PseudoTensor& T1, const PseudoTensor& T2);

      bool is_overlap(const PseudoTensor& T1, const PseudoTensor& T2);

      PseudoTensor psesudocontract(const PseudoTensor& T1, const PseudoTensor& T2);

      float get_cost(const PseudoTensor& T1, const PseudoTensor& T2);

      ary2d_pten tensor_set;

      float xi_min;

      int numT;

      char *netorder;

  };

};
#endif
