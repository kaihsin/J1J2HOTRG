#include "uni10/uni10_type.h"
#include "uni10/uni10_api/network_tools/network_tools.h"

namespace uni10{

  std::vector<_Swap> recSwap(std::vector<uni10_int32>& _ord) { //Given the reshape order out to in.
    //int ordF[n];
    uni10_int32 n = _ord.size();
    std::vector<uni10_int32> ordF(n);
    for(uni10_int32 i = 0; i < n; i++)
      ordF[i] = i;
    return recSwap(_ord, ordF);
  }

  std::vector<_Swap> recSwap(std::vector<uni10_int32>& _ord, std::vector<uni10_int32>& ordF) { //Given the reshape order out to in.

    uni10_int32 n = _ord.size();
    std::vector<uni10_int32> ord = _ord;
    std::vector<_Swap> swaps;
    _Swap sg;
    uni10_int32 tmp;
    for(uni10_int32 i = 0; i < n - 1; i++)
      for(uni10_int32 j = 0; j < n - i - 1; j++)
        if(ord[j] > ord[j + 1]) {
          sg.b1 = ordF[ord[j + 1]];
          sg.b2 = ordF[ord[j]];
          tmp = ord[j];
          ord[j] = ord[j + 1];
          ord[j + 1] = tmp;
          swaps.push_back(sg);
        }

    return swaps;

  }

};
