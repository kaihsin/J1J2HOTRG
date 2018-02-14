#ifndef __UNI10_NETTOOLS_H__
#define __UNI10_NETTOOLS_H__

#include <vector>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

namespace uni10{

  typedef struct{
    int b1;
    int b2;
  }_Swap;

  std::vector<_Swap> recSwap(std::vector<int>& ord, std::vector<int>& ordF);

  std::vector<_Swap> recSwap(std::vector<int>& ord);	//Given the reshape order out to in.

  // trim from start
  static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
  }
  // trim from end
  static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
  }
  // trim from both ends
  static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
  }

}

#endif
