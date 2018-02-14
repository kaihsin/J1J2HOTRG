/****************************************************************************
*  @file uni10/uni10_api/node.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    National Taiwan University
*    National Tsing-Hua University

*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief Generic header file for uni10 data strucutres
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef __UNI10_NODE_DEV_H__
#define __UNI10_NODE_DEV_H__

#include <utility>

#include "uni10/uni10_api/UniTensor.h"
#include "uni10/uni10_api/hirnk_linalg_inplace.h"
#include "uni10/uni10_api/network_tools/network_tools.h"

namespace uni10{

  class Bond;

  class Node_dev;

  std::ostream& operator<< (std::ostream& os, const Node_dev& nd);

  class Node_dev {

    public:

      Node_dev();

      ~Node_dev();

      void init_swap_arr();

      void init_Tout();

      void prepare4cout( std::vector<std::vector<std::vector<std::string> > >& info4cout, std::vector<uni10_uint64>& ln_pos, 
          std::map<std::string, std::map<uni10_int, uni10_uint64> >& other_nd_lableDim ) const;

      void get_cout_strs(std::vector<std::string>& out_strs, const std::vector<std::vector<std::vector<std::string> > >& info4cout, const std::vector<uni10_uint64>& ln_pos ) const;

      void merge();

      std::vector<uni10_uint64> cal_elemNums(std::map<std::string, std::map<uni10_int, uni10_uint64> >& other_nd_lableDim) const;

      friend std::ostream& operator<< (std::ostream& os, const Node_dev& nd);

      friend class Network_dev;

      friend class layer;

    private:

      std::string nd_name;

      std::vector< std::string* > tens_names;

      std::vector< std::vector<uni10_int>* > tens_labels;

      // void pointer of a tensor and its type.
      // The index of tensors in the vector is following the contraction order.
      std::vector<std::pair<const void*, int> > tens_in_nd;

      UniTensor<uni10_double64>* Tr;

      UniTensor<uni10_complex128>* Tc;

      std::pair<void*, int> Tout;

      // Prepare for cout
      std::vector<uni10_uint64> elemNum_per_step;

      // wether the tensors in the node are the same type or not.
      // Recording the length of the continuous real tensors in the list if it is mixed.
      std::pair<bool, uni10_int> isMix_TrLen;

      // addSwap
      std::vector< uni10_bool* > swapflags;

      std::vector< std::vector<_Swap>* > swaps_arr;

  };

}; /* namespace uni10 */

#endif
