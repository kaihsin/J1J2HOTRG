/****************************************************************************
 *  @file Network.cpp
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
 *  @brief Implementation file for Node and Network classes 
 *  @author Yun-Hsuan Chou, Ying-Jer Kao
 *  @date 2017-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/
#include <algorithm>

#include "uni10/uni10_api/network_tools_dev/layer.h"

namespace uni10{

    std::ostream& operator<< (std::ostream& os, const layer& ly){

      std::vector< std::vector<std::string> > each_nd_out_strs(ly.nds_in_layer.size());
      std::map<std::string, std::map<uni10_int, uni10_uint64> > other_nd_labelDim;
      std::vector<std::vector<std::vector<std::string> > > info4cout;
      std::vector<uni10_uint64> ln_pos;

      for(uni10_uint64 l = 0; l < ly.nds_in_layer.size(); l++){
        
        ly.nds_in_layer[l]->prepare4cout(info4cout, ln_pos, other_nd_labelDim);
        ly.nds_in_layer[l]->get_cout_strs(each_nd_out_strs[l], info4cout, ln_pos);
        info4cout.clear();
        ln_pos.clear();

      }

      // Combine
      uni10_uint64 max_line = each_nd_out_strs[0].size();

      for(uni10_uint64 k = 0; k < each_nd_out_strs.size()-1; k++)
        max_line = (max_line < each_nd_out_strs[k+1].size()) ? each_nd_out_strs[k+1].size() : max_line;

      for(uni10_uint64 n = 0; n < each_nd_out_strs.size(); n++){

        for(uni10_uint64 m = 0; m < max_line - each_nd_out_strs[n].size(); m++){

          each_nd_out_strs[n].push_back(std::string(each_nd_out_strs[n][0].size(), ' '));

        }

      }

      std::vector<std::string> out_strs(max_line);

      for(uni10_uint64 c = 0; c < max_line; c++){
        out_strs[c] = each_nd_out_strs[0][c];
        for(uni10_uint64 d = 0; d < each_nd_out_strs.size()-1; d++)
          out_strs[c] += each_nd_out_strs[d+1][c];

      }

      fprintf(stdout, "\n");
      for(uni10_uint64 j = 0; j < out_strs.size(); j++)
        fprintf(stdout, "%s\n", out_strs[j].c_str());
      fprintf(stdout, "\n");

      os << "";

      return os;

    }

    layer::layer(){

    }

    layer::~layer(){

    }

    void layer::merge(){

      uni10_int i;

      for(i = 0; i < (uni10_int)nds_in_layer.size(); i++){

        nds_in_layer[i]->merge();

      }

    }

    void layer::get_cout_strs(std::vector<std::string>& out_strs, std::map<std::string, std::map<uni10_int, uni10_uint64> >& other_nd_labelDim) const{

      std::vector< std::vector<std::string> > each_nd_out_strs(nds_in_layer.size());
      
      std::vector<std::vector<std::vector<std::string> > > info4cout;
      std::vector<uni10_uint64> ln_pos;

      for(uni10_uint64 l = 0; l < nds_in_layer.size(); l++){
        
        nds_in_layer[l]->prepare4cout(info4cout, ln_pos, other_nd_labelDim);
        nds_in_layer[l]->get_cout_strs(each_nd_out_strs[l], info4cout, ln_pos);
        info4cout.clear();
        ln_pos.clear();

      }

      // Combine
      uni10_uint64 max_line = each_nd_out_strs[0].size();

      for(uni10_uint64 k = 0; k < each_nd_out_strs.size()-1; k++)
        max_line = (max_line < each_nd_out_strs[k+1].size()) ? each_nd_out_strs[k+1].size() : max_line;

      for(uni10_uint64 n = 0; n < each_nd_out_strs.size(); n++){

        uni10_uint64 ori_size = each_nd_out_strs[n].size();

        for(uni10_uint64 m = 0; m < max_line - ori_size; m++){

          each_nd_out_strs[n].push_back(std::string(each_nd_out_strs[n][0].size(), ' '));

        }

      }

      out_strs.assign(max_line, std::string());

      for(uni10_uint64 c = 0; c < max_line; c++){
        out_strs[c] = each_nd_out_strs[0][c];
        for(uni10_uint64 d = 0; d < each_nd_out_strs.size()-1; d++)
          out_strs[c] += each_nd_out_strs[d+1][c];

      }

    }

}; /* namespace uni10 */
