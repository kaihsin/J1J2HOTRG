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
#ifndef __UNI10_LAYER_H__
#define __UNI10_LAYER_H__

#include "uni10/uni10_api/UniTensor.h"
#include "uni10/uni10_api/network_tools_dev/node.h"

namespace uni10{

  class Node_dev;

  class layer;

  std::ostream& operator<< (std::ostream& os, const layer& ly);

  class layer {

    public:

      layer();

      ~layer();

      void merge();

      void get_cout_strs(std::vector<std::string>& out_strs, std::map<std::string, std::map<uni10_int, uni10_uint64> >& other_nd_labelDim) const;

      friend std::ostream& operator<< (std::ostream& os, const layer& ly);

      friend class Network_dev;

      friend class Node_dev;

    private:

      const std::vector<std::string> *names;

      std::vector<Node_dev*> nds_in_layer;

  };

}; /* namespace uni10 */

#endif
