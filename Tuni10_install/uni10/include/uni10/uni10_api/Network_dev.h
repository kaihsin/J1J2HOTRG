/****************************************************************************
 *  @file Network.h
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
 *  @brief Header file for Newtork class
 *  @author Yun-Da Hsieh
 *  @date 2014-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/
#ifndef NETWORK_DEV_H
#define NETWORK_DEV_H

#include <utility>

#include "uni10/uni10_api/UniTensor.h"
#include "uni10/uni10_api/hirnk_linalg_inplace.h"
#include "uni10/uni10_api/network_tools_dev/layer.h"
#include "uni10/uni10_api/network_tools_dev/NetOrder.h"

namespace uni10 {

  template<typename T>
    class UniTensor;

  class layer;

  class Node_dev;

  class Network_dev;

  std::ostream& operator<< (std::ostream& os, const Network_dev& net);

  class Network_dev {

    public:

      Network_dev(const std::string& fname);

      ~Network_dev();

      void pre_construct(bool force = true);

      void putTensor(int idx, const UniTensor<uni10_double64>& UniT, bool force = true);

      void putTensor(int idx, const UniTensor<uni10_complex128>& UniT, bool force = true);

      void putTensor(const std::string& name, const UniTensor<uni10_double64>& UniT, bool force = true);

      void putTensor(const std::string& name, const UniTensor<uni10_complex128>& UniT, bool force = true);

      void putTensorT(int idx, const UniTensor<uni10_double64>& UniT, bool force = true);

      void putTensorT(int idx, const UniTensor<uni10_complex128>& UniT, bool force = true);

      void putTensorT(const std::string& name, const UniTensor<uni10_double64>& UniT, bool force = true);

      void putTensorT(const std::string& name, const UniTensor<uni10_complex128>& UniT, bool force = true);

      void putTensorD(int idx, const UniTensor<uni10_double64>& UniT, bool force = true);

      void putTensorD(int idx, const UniTensor<uni10_complex128>& UniT, bool force = true);

      void putTensorD(const std::string& name, const UniTensor<uni10_double64>& UniT, bool force = true);

      void putTensorD(const std::string& name, const UniTensor<uni10_complex128>& UniT, bool force = true);

      void launch(UniTensor<uni10_double64>& Tout, const std::string& Tname="");

      void launch(UniTensor<uni10_complex128>& Tout, const std::string& Tname="");

      std::string get_contract_order() const;

      friend std::ostream& operator<< (std::ostream& os, const Network_dev& net);

      friend class layer;

      friend class Node_dev;

      template<typename T>
        friend class UniTensor;

    private:

      // Name of the network file which can be written manually or generated from UNI10 GUI.
      std::string fname;

      //whether or not the network is ready for contraction, construct=> load=true, destruct=>load=false
      bool load, isChanged; 

      // Tensor List.
      // For saving the pointers of UniTensors which are put in.
      // Under normal circumstance, these pointers have to keep in constant.
      std::vector< std::pair<const void*, int> > tensor_type;

      // whether or not the network has specific order.
      bool ordered; 

      // Contraction order.
      std::vector<uni10_int> contract_order;

      // labels correspond to the tensor in the list.
      std::vector< std::vector<uni10_int> > label_arr;

      std::vector< uni10_int > iBondNums;

      // The name of each tenosr.
      std::vector<std::string> names;

      // To count the position index of each tensor in the network file.
      std::map<std::string, std::vector<uni10_int> > name2pos;

      void fromfile(const std::string& fname);

      // The structure of contraction B-tree.
      std::vector<layer*> contraction_layers;

      // pos2pos_int_tre[ten_idx_in_tensor_list][0] = layer idx, l.
      // pos2pos_int_tre[ten_idx_in_tensor_list][1] = node idx, n, in layer l.
      // pos2pos_int_tre[ten_idx_in_tensor_list][2] = idx, i, which is the idx of contraction order in node, n,.
      std::vector< std::vector<uni10_int> > pos2pos_in_tree;

      void putTensor_driver(int idx, const void* UniT, int typeID, bool force = true);

      // Initialize each layers in B-tree, which is a recursive function.
      void init_layers(uni10_int& idx_ly, std::vector<uni10_int>& _contrac_order, std::vector<uni10_int>& ori_ten_idx,
          std::vector<std::pair<const void*, int> >& _tens_type, std::vector<std::vector<uni10_int>* >& _label_arr, std::vector<std::string* >& _namas);
   
      // Construct contraction B-tree.
      void construct();

      // Destruct contraction B-tree.
      void destruct();

      void destruct_tree_only();

      void init_swap_arr_in_nodes();

      void del_swap_arr_in_nodes();

      void addSwap();

  };

};  /* namespace uni10 */

#endif /* NETWORK<uni10_type>_H */
