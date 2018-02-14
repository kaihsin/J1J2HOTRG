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

#include "uni10/uni10_error.h"
#include "uni10/uni10_api/network_tools_dev/node.h"

namespace uni10{

    std::ostream& operator<< (std::ostream& os, const Node_dev& nd){

      std::map<std::string, std::map<uni10_int, uni10_uint64> > other_nd_labelDim;
      std::vector<std::vector<std::vector<std::string> > > info4cout;
      std::vector<uni10_uint64> ln_pos;

      nd.prepare4cout(info4cout, ln_pos, other_nd_labelDim);
     
      std::vector<std::string> out_strs;
      nd.get_cout_strs(out_strs, info4cout, ln_pos);

      fprintf(stdout, "\n");
      for(uni10_uint64 l = 0; l < out_strs.size(); l++)
        fprintf(stdout, "%s\n", out_strs[l].c_str());
      fprintf(stdout, "\n");

      os << "";

      return os;

    }

    Node_dev::Node_dev(): nd_name(""), Tr(NULL), Tc(NULL), isMix_TrLen(std::pair<bool, int>(false, 0)){

    }

    Node_dev::~Node_dev(){

      if(Tr != NULL){

        delete [] Tr;

      }
      else if (Tc != NULL){

        delete [] Tc;

      }

    }

    void Node_dev::init_swap_arr(){

      uni10_error_msg(!Qnum::isFermionic(), "%s", "Unexpected error. Please connect to the developers of uni10.");

    }

    void Node_dev::init_Tout(){

      uni10_error_msg(Tr != NULL && Tc != NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );

      uni10_uint64 i;

      bool isTc = false;

      uni10_int cnt_Tr = 0;

      for( i = 0; i < tens_in_nd.size(); i++){

        if(tens_in_nd[i].second == 2)
          isTc = true;

        if(tens_in_nd[i].second == 1 && !isTc)
          cnt_Tr++;
       
        if(i == 0)
          continue;

        else if(tens_in_nd[i-1].second != tens_in_nd[i].second){

          isMix_TrLen = std::pair<bool, uni10_int>(true, cnt_Tr);

        }

      }

      if(!isTc){

        Tr = new UniTensor<uni10_double64>[1];
        Tout = std::pair<const void*, int>((void*)Tr, Tr->typeID());

      }
      else{

        Tc = new UniTensor<uni10_complex128>[1];
        Tout = std::pair<const void*, int>((void*)Tc, Tc->typeID());

      }

    }

    void Node_dev::prepare4cout(std::vector<std::vector<std::vector<std::string> > >& info4cout, std::vector<uni10_uint64>& ln_pos
        , std::map<std::string, std::map<uni10_int, uni10_uint64> >& other_nd_labelDim ) const{

      std::vector<std::string>   tens_name_str;
      std::vector<std::string>   tens_elemNum_str;

      std::vector<uni10_uint64>  _ln_pos(tens_in_nd.size());

      std::vector<std::string> nd_info(2, std::string());
      std::vector<std::vector<std::string> > sub_info;

      //
      // space_len[0]: The number of space before the names of tensor. 
      // space_len[1]: The number of space before the element number of tensor. 
      //
      std::vector<uni10_uint64>  space_len(2, 0);

      std::vector<uni10_uint64> leaf_tens_name_pos;
      std::vector<uni10_uint64> leaf_tens_elemNum_pos;
      
      std::vector<uni10_uint64> conDims = this->cal_elemNums(other_nd_labelDim);
      std::vector<std::string> conNames(tens_in_nd.size()-1);
      std::vector<std::string> conDims_str;

      uni10_uint64 i;

      for(i = 0; i < conNames.size(); i++ ){
        conNames[conNames.size()-1-i] = (i != conNames.size()-1) ? *tens_names[i+1]+"\'" : nd_name;
        conDims_str.push_back("(" + std::to_string(conDims[i]) +")");
      }

      // Fill spaces
      //
      uni10_uint64 longer;
      
      for(i = 0; i < conNames.size(); i++){

        longer= std::max(conNames[i].size(), conDims_str[i].size());

        space_len[0] = conDims_str[i].size() <= conNames[i].size() ? 0 : (conDims_str[i].size() - conNames[i].size())/2;
        space_len[1] = conDims_str[i].size() <  conNames[i].size() ? (conNames[i].size() - conDims_str[i].size())/2 : 0;

        conNames[i]    = std::string(space_len[0], ' ') + conNames[i]    + std::string(longer-space_len[0]-conNames[i].size(), ' ');
        conDims_str[i] = std::string(space_len[1], ' ') + conDims_str[i] + std::string(longer-space_len[1]-conDims_str[i].size(), ' ') ;

      }

      for(i = 0; i < conNames.size(); i++ ){
        nd_info[0] = " " + conNames[i] + " ";
        nd_info[1] = " " + conDims_str[i] + " ";
        sub_info.push_back(nd_info);
        info4cout.push_back(sub_info);
        sub_info.clear();
      }

      std::vector<uni10_uint64> name_offset;

      uni10_uint64 ten_elemNum;

      for(i = 0; i < tens_in_nd.size(); i++){

        if(tens_labels[i] != NULL){

          if(tens_in_nd[i].second == 1) 
            ten_elemNum = ((UniTensor<uni10_double64>*)tens_in_nd[i].first)->elemNum();
          else if(tens_in_nd[i].second == 2)
            ten_elemNum = ((UniTensor<uni10_complex128>*)tens_in_nd[i].first)->elemNum();

        }else{

          std::map<uni10_int, uni10_uint64> tmp_labelDim = other_nd_labelDim[*tens_names[i]];
          ten_elemNum = 1;
          std::map<uni10_int, uni10_uint64>::const_iterator it = tmp_labelDim.begin();
          for(; it != tmp_labelDim.end(); ++it)
            ten_elemNum *= it->second;

        }

        tens_name_str.push_back(*tens_names[i]);
        tens_elemNum_str.push_back("(" + std::to_string(ten_elemNum) + ")");

        space_len[0] = tens_elemNum_str.back().size() <= tens_name_str.back().size() ? 0 : (tens_elemNum_str.back().size() - tens_name_str.back().size())/2;
        space_len[1] = tens_elemNum_str.back().size() <  tens_name_str.back().size() ? (tens_name_str.back().size() - tens_elemNum_str.back().size())/2 : 0;

        longer= std::max(tens_name_str[i].size(), tens_elemNum_str[i].size());

        tens_name_str.back()    = std::string(space_len[0]+1, ' ') + tens_name_str.back() + std::string(longer-space_len[0]-tens_name_str.back().size()+1, ' ');
        tens_elemNum_str.back() = std::string(space_len[1]+1, ' ') + tens_elemNum_str.back() + std::string(longer-space_len[1]-tens_elemNum_str.back().size()+1, ' ') ;
 
        nd_info[0] = tens_name_str.back();
        nd_info[1] = tens_elemNum_str.back();
        sub_info.push_back(nd_info);

      }

      info4cout.push_back(sub_info);
      sub_info.clear();

      // Find _ln_pos[0]
      uni10_uint64 max_name_str = info4cout[0][0][0].size();
      uni10_uint64 max_elemNum_str = info4cout[0][0][1].size();

      for(uni10_uint64 z = 0; z < info4cout.size()-1; z++){

        max_name_str = (max_name_str < info4cout[z+1][0][0].size()) ? info4cout[z+1][0][0].size() : max_name_str;
        max_elemNum_str = (max_elemNum_str < info4cout[z+1][0][1].size()) ? info4cout[z+1][0][1].size() : max_elemNum_str;

      }

      uni10_uint64 max_str = max_name_str > max_elemNum_str ? max_name_str : max_elemNum_str;
      _ln_pos[0] = max_str % 2 == 0 ? max_str/2 - 1 : max_str/2;

      // Filled space.
      for(uni10_uint64 z = 0; z < info4cout.size(); z++){

        //uni10_uint64 front_name = (max_str-info4cout[z][0][0].size()) == 1 ?  1 : (max_str-info4cout[z][0][0].size())/2;
        uni10_uint64 front_name = (max_str-info4cout[z][0][0].size())/2;
        uni10_uint64 end_name   = max_str-front_name-info4cout[z][0][0].size();

        //uni10_uint64 front_elemNum = (max_str-info4cout[z][0][1].size()) == 1 ? 1 : (max_str-info4cout[z][0][1].size())/2;
        uni10_uint64 front_elemNum = (max_str-info4cout[z][0][1].size())/2;
        uni10_uint64 end_elemNum   = max_str-front_elemNum-info4cout[z][0][1].size();

        info4cout[z][0][0] = std::string(front_name, ' ') + info4cout[z][0][0] + std::string(end_name, ' ');
        info4cout[z][0][1] = std::string(front_elemNum, ' ') + info4cout[z][0][1] + std::string(end_elemNum, ' ');
      }

      //Find the rest of _ln_pos.
      uni10_uint64 offset = info4cout.back()[0][0].size();
      uni10_uint64 mid;

      for(uni10_uint64 z = 0; z < info4cout.back().size()-1; z++){

        mid = info4cout.back()[z+1][0].size() % 2 == 0 ? (info4cout.back()[z+1][0].size()/2) - 1: info4cout.back()[z+1][0].size()/2;
        _ln_pos[z+1] = offset + mid;
        offset += info4cout.back()[z+1][0].size();

      }
     
      ln_pos = _ln_pos;

    }

    void Node_dev::get_cout_strs(std::vector<std::string>& out_strs, const std::vector<std::vector<std::vector<std::string> > >& info4cout, const std::vector<uni10_uint64>& ln_pos) const{

      uni10_uint64 line_len = 0;
      for(uni10_uint64 l = 0; l < info4cout.back().size(); l++)
        line_len += info4cout.back()[l][0].size();

      uni10_uint64 line_num = ( 2*info4cout.size()-1 ) * 2;
      out_strs.assign(line_num, std::string(line_len, ' '));

      for(uni10_uint64 i = 0; i < line_num; i++){

        uni10_uint64 ly_in_nd_idx = i/4;

        if(i % 4 == 0){

          if(info4cout[ly_in_nd_idx].size() != 1){
            out_strs[i] = info4cout[ly_in_nd_idx][0][0];

            for(uni10_uint64 k = 0; k < info4cout[ly_in_nd_idx].size()-1; k++)
              out_strs[i] += info4cout[ly_in_nd_idx][k+1][0];

          }else{

            for(uni10_uint64 k = 0; k < line_len; k++){

              if(k < info4cout[ly_in_nd_idx][0][0].size() - 1){

                out_strs[i][k] = info4cout[ly_in_nd_idx][0][0][k];

              }else{

                if(ly_in_nd_idx != 0){

                  for(uni10_uint64 t = 0; t < ly_in_nd_idx; t++)
                    out_strs[i][ln_pos[ln_pos.size() - 1 - t]] = '|';

                }

              }

            }

          }

        }

        else if(i % 4 == 1){

          if(info4cout[ly_in_nd_idx].size() != 1){
            out_strs[i] = info4cout[ly_in_nd_idx][0][1];

            for(uni10_uint64 k = 0; k < info4cout[ly_in_nd_idx].size()-1; k++)
              out_strs[i] += info4cout[ly_in_nd_idx][k+1][1];
          }else{

            for(uni10_uint64 k = 0; k < line_len; k++){

              if(k < info4cout[ly_in_nd_idx][0][1].size() - 1){

                out_strs[i][k] = info4cout[ly_in_nd_idx][0][1][k];

              }else{

                if(k <= ln_pos[ln_pos.size() - 1 - ly_in_nd_idx])
                  out_strs[i][k] = '_';

                if(ly_in_nd_idx != 0 ){

                  for(uni10_uint64 t = 0; t < ly_in_nd_idx; t++)
                    out_strs[i][ln_pos[ln_pos.size() - 1 - t]] = '|';
                  
                }

              }

            }

          }

        }

        else{

          out_strs[i][ln_pos[0]] = '|';
          for(uni10_uint64 t = 0; t < ly_in_nd_idx+1; t++)
            out_strs[i][ln_pos[ln_pos.size() - 1 - t]] = '|';

        }

      }

    }

    std::vector<uni10_uint64> Node_dev::cal_elemNums(std::map<std::string, std::map<uni10_int, uni10_uint64> >& other_nd_labelDim) const{

      std::map<uni10_int, uni10_uint64> label_dim;
      std::map<std::string, std::vector<uni10_int> > other_nd_label;

      uni10_uint64 n;

      for(n = 0; n < tens_in_nd.size(); n++){

        uni10_uint64 b;

        if(tens_labels[n] != NULL){

          for(b = 0; b < tens_labels[n]->size(); b++){

            if(tens_in_nd[n].second == 1)
              label_dim[(*tens_labels[n])[b]] = ((UniTensor<uni10_double64>*)tens_in_nd[n].first)->bond(b).dim();
            else
              label_dim[(*tens_labels[n])[b]] = ((UniTensor<uni10_double64>*)tens_in_nd[n].first)->bond(b).dim();

          }

        }else{
          
          std::map<uni10_int, uni10_uint64> _label_dim = other_nd_labelDim[*tens_names[n]];

          std::map<uni10_int, uni10_uint64>::const_iterator it = _label_dim.begin();

          for(; it != _label_dim.end(); it++){
            label_dim[it->first] = it->second;
            other_nd_label[*tens_names[n]].push_back(it->first);

          }

        }

      }

      // Find matched labels.
      std::vector<uni10_uint64> dims(tens_in_nd.size()-1);
      
      std::vector<uni10_int> remained_label = tens_labels[0] != NULL ? *tens_labels[0] : other_nd_label[*tens_names[0]];
      std::vector<uni10_int> match_label, tmp_label;


      std::vector<uni10_int>::iterator it;

      for(n = 0; n < tens_in_nd.size()-1; n++){

        std::vector<uni10_int> tmp_next_ten = tens_labels[n+1] != NULL ? *tens_labels[n+1] : other_nd_label[*tens_names[n+1]];

        uni10_uint64 m;

        for( m = 0; m < remained_label.size(); m++){

          it = std::find(tmp_next_ten.begin(), tmp_next_ten.end(), remained_label[m]);

          if(it == tmp_next_ten.end())
            tmp_label.push_back(remained_label[m]);
          else
            match_label.push_back(remained_label[m]);

        }

        remained_label = tmp_label;

        uni10_uint64 k; 


        for(k = 0; k < tmp_next_ten.size(); k++){

          it = std::find(match_label.begin(), match_label.end(), tmp_next_ten[k]);

          if(it == match_label.end())
            remained_label.push_back(tmp_next_ten[k]);

        }

        if(n == tens_in_nd.size()-2){

          std::map<uni10_int, uni10_uint64> Tout_labelDim;

          for(uni10_uint64 v = 0; v < remained_label.size(); v++)
            Tout_labelDim[remained_label[v]] = label_dim[remained_label[v]];

          other_nd_labelDim[nd_name] = Tout_labelDim;

        }
        
        uni10_uint64 dim = 1;

        for(m = 0; m < remained_label.size(); m++)
          dim *= label_dim[remained_label[m]];

        dims[dims.size()-1-n] = dim;

        tmp_label.clear();

      }

      return dims;
    
    }

    void Node_dev::merge(){

      uni10_error_msg(Qnum::isFermionic() ^ swapflags.size(), "%s", "Unexpected error. Please connect to the developers of uni10.");
      uni10_error_msg(Qnum::isFermionic() ^ swaps_arr.size(), "%s", "Unexpected error. Please connect to the developers of uni10.");

      if(isMix_TrLen.first){
        
        UniTensor<uni10_complex128> tmpT1(*(UniTensor<uni10_complex128>*)tens_in_nd[isMix_TrLen.second].first);

        uni10_uint64 complex1st = isMix_TrLen.second;

        if(tens_labels[complex1st] == NULL){

          ((UniTensor<uni10_complex128>*)tens_in_nd[complex1st].first)->clear();

        }
        else{

          tmpT1.setLabel(*tens_labels[complex1st]);
          // AddSwaps
          uni10_error_msg(Qnum::isFermionic() && swapflags[complex1st] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
          uni10_error_msg(Qnum::isFermionic() && swaps_arr[complex1st] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
          if(Qnum::isFermionic() && !swapflags[complex1st]){

            tmpT1.addGate(*swaps_arr[complex1st]);
            *swapflags[complex1st] = true;

          }

        }

        if(isMix_TrLen.second > 1){

          UniTensor<uni10_double64> Tr1(*(UniTensor<uni10_double64>*)tens_in_nd[0].first);
          UniTensor<uni10_double64> Tr2(*(UniTensor<uni10_double64>*)tens_in_nd[1].first);
          UniTensor<uni10_double64> Tr3;

          if(tens_labels[0] == NULL){

            ((UniTensor<uni10_double64>*)tens_in_nd[0].first)->clear();

          }
          else{

            Tr1.setLabel(*tens_labels[0]);
            // AddSwaps
            uni10_error_msg(Qnum::isFermionic() && swapflags[0] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            uni10_error_msg(Qnum::isFermionic() && swaps_arr[0] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            if(Qnum::isFermionic() && !swapflags[0]){

              Tr1.addGate(*swaps_arr[0]);
              *swapflags[0] = true;

            }

          }

          if(tens_labels[1] == NULL){

            ((UniTensor<uni10_double64>*)tens_in_nd[1].first)->clear();

          }
          else{

            Tr2.setLabel(*tens_labels[1]);
            // AddSwaps
            uni10_error_msg(Qnum::isFermionic() && swapflags[1] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            uni10_error_msg(Qnum::isFermionic() && swaps_arr[1] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            if(Qnum::isFermionic() && !swapflags[1]){

              Tr2.addGate(*swaps_arr[1]);
              *swapflags[1] = true;

            }

          }

          contract_driver[0]((void*)&Tr3, (void*)&Tr1, (void*)&Tr2, true);

          uni10_int i = 2;

          for( i = 2; i < isMix_TrLen.second; i++){

            Tr2 = *(UniTensor<uni10_double64>*)tens_in_nd[i].first;

            if(tens_labels[i] == NULL){

              ((UniTensor<uni10_double64>*)tens_in_nd[i].first)->clear();

            }
            else{

              Tr2.setLabel(*tens_labels[i]);
              // AddSwaps
              uni10_error_msg(Qnum::isFermionic() && swapflags[i] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              uni10_error_msg(Qnum::isFermionic() && swaps_arr[i] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              if(Qnum::isFermionic() && !swapflags[i]){

                Tr2.addGate(*swaps_arr[i]);
                *swapflags[i] = true;

              }

            }

            if(i % 2 == 0)
              contract_driver[0]((void*)&Tr1, (void*)&Tr3, (void*)&Tr2, true);
            else
              contract_driver[0]((void*)&Tr3, (void*)&Tr1, (void*)&Tr2, true);

          }

          if(i % 2 == 0){
            contract_driver[2](Tout.first, (void*)&Tr3, (void*)&tmpT1, true);
          }

          if(i % 2 == 1){
            contract_driver[2](Tout.first, (void*)&Tr1, (void*)&tmpT1, true);
          }

        }else{

          uni10_int rest_idx = isMix_TrLen.second == 0 ? 1 : 0;

          if(tens_in_nd[rest_idx].second == 1){

            UniTensor<uni10_double64> tmp_rest(*(UniTensor<uni10_double64>*)tens_in_nd[rest_idx].first);

            if(tens_labels[rest_idx] == NULL){

              ((UniTensor<uni10_double64>*)tens_in_nd[rest_idx].first)->clear();

            }

            else{

              tmp_rest.setLabel(*tens_labels[rest_idx]);
              // AddSwaps
              uni10_error_msg(Qnum::isFermionic() && swapflags[rest_idx] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              uni10_error_msg(Qnum::isFermionic() && swaps_arr[rest_idx] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              if(Qnum::isFermionic() && !swapflags[rest_idx]){

                tmp_rest.addGate(*swaps_arr[rest_idx]);
                *swapflags[rest_idx] = true;

              }

            }

            contract_driver[1](Tout.first, (void*)&tmpT1, (void*)&tmp_rest, true);

          }else{

            UniTensor<uni10_complex128> tmp_rest(*(UniTensor<uni10_complex128>*)tens_in_nd[rest_idx].first);


            if(tens_labels[rest_idx] == NULL){

              ((UniTensor<uni10_complex128>*)tens_in_nd[rest_idx].first)->clear();

            }
            else{

              tmp_rest.setLabel(*tens_labels[rest_idx]);
              // AddSwaps
              uni10_error_msg(Qnum::isFermionic() && swapflags[rest_idx] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              uni10_error_msg(Qnum::isFermionic() && swaps_arr[rest_idx] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              if(Qnum::isFermionic() && !swapflags[rest_idx]){

                tmp_rest.addGate(*swaps_arr[rest_idx]);
                *swapflags[rest_idx] = true;

              }

            }

            contract_driver[3](Tout.first, (void*)&tmpT1, (void*)&tmp_rest, true);

          }

        }

        uni10_uint64 offset = (isMix_TrLen.second > 1) ? isMix_TrLen.second + 1 : 2;
        uni10_uint64 j = 0;

        for( ; j < tens_in_nd.size()-offset; j++){

          if(tens_in_nd[j+offset].second == 1){

            UniTensor<uni10_double64> tmpT2(*(UniTensor<uni10_double64>*)tens_in_nd[j+offset].first);

            if(tens_labels[j+offset] == NULL){

              ((UniTensor<uni10_double64>*)tens_in_nd[j+offset].first)->clear();

            }
            else{

              tmpT2.setLabel(*tens_labels[j+offset]);
              // AddSwaps
              uni10_error_msg(Qnum::isFermionic() && swapflags[j+offset] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              uni10_error_msg(Qnum::isFermionic() && swaps_arr[j+offset] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              if(Qnum::isFermionic() && !swapflags[j+offset]){

                tmpT2.addGate(*swaps_arr[j+offset]);
                *swapflags[j+offset] = true;

              }

            }

            if(j % 2 == 0)
              contract_driver[1]((void*)&tmpT1, Tout.first, (void*)&tmpT2, true);
            else
              contract_driver[1](Tout.first, (void*)&tmpT1, (void*)&tmpT2, true);

          }

          else if (tens_in_nd[j+offset].second == 2){

            UniTensor<uni10_complex128> tmpT2(*(UniTensor<uni10_complex128>*)tens_in_nd[j+offset].first);

            if(tens_labels[j+offset] == NULL){

              ((UniTensor<uni10_complex128>*)tens_in_nd[j+offset].first)->clear();

            }
            else{

              tmpT2.setLabel(*tens_labels[j+offset]);
              // AddSwaps
              uni10_error_msg(Qnum::isFermionic() && swapflags[j+offset] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              uni10_error_msg(Qnum::isFermionic() && swaps_arr[j+offset] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              if(Qnum::isFermionic() && !swapflags[j+offset]){

                tmpT2.addGate(*swaps_arr[j+offset]);
                *swapflags[j+offset] = true;

              }

            }

            if(j % 2 == 0)
              contract_driver[3]((void*)&tmpT1, Tout.first, (void*)&tmpT2, true);
            else
              contract_driver[3](Tout.first, (void*)&tmpT1, (void*)&tmpT2, true);

          }else{

            uni10_error_msg(true, "%s", "Unexpected error. Please connect to the developers of uni10.");

          }

        }

        if(j % 2 == 1){

          *(UniTensor<uni10_complex128>*)Tout.first = tmpT1;

        }

      }

      else{

        // Tensors in contraction list have same type.

        int driver_type;

        if(Tout.second == 1){

          driver_type = 0;

          UniTensor<uni10_double64> tmpT1(*(UniTensor<uni10_double64>*)tens_in_nd[0].first);
          UniTensor<uni10_double64> tmpT2(*(UniTensor<uni10_double64>*)tens_in_nd[1].first);

          // Set labels and add Swap gate into T1.
          if(tens_labels[0] == NULL){

            ((UniTensor<uni10_double64>*)tens_in_nd[0].first)->clear();

          }
          else{

            tmpT1.setLabel(*tens_labels[0]);
            // AddSwaps
            uni10_error_msg(Qnum::isFermionic() && swapflags[0] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            uni10_error_msg(Qnum::isFermionic() && swaps_arr[0] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            if(Qnum::isFermionic() && !swapflags[0]){

             tmpT1.addGate(*swaps_arr[0]);
             *swapflags[0] = true;

            }

          }

          // Set labels and add Swap gate into T2.
          
          if(tens_labels[1] == NULL){

            ((UniTensor<uni10_double64>*)tens_in_nd[1].first)->clear();

          }
          else{

            tmpT2.setLabel(*tens_labels[1]);
            // AddSwaps
            uni10_error_msg(Qnum::isFermionic() && swapflags[1] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            uni10_error_msg(Qnum::isFermionic() && swaps_arr[1] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            if(Qnum::isFermionic() && !swapflags[1]){

             tmpT2.addGate(*swaps_arr[1]);
             *swapflags[1] = true;

            }

          }

          contract_driver[driver_type](Tout.first, (void*)&tmpT1, (void*)&tmpT2, true);

          uni10_uint64 i;

          for( i = 2; i < tens_in_nd.size(); i++){

            tmpT2 = *(UniTensor<uni10_double64>*)tens_in_nd[i].first;

            if(tens_labels[i] == NULL){

              ((UniTensor<uni10_double64>*)tens_in_nd[i].first)->clear();

            }
            else{

              tmpT2.setLabel(*tens_labels[i]);
              // AddSwaps
              uni10_error_msg(Qnum::isFermionic() && swapflags[i] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              uni10_error_msg(Qnum::isFermionic() && swaps_arr[i] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              if(Qnum::isFermionic() && !swapflags[i]){

                tmpT2.addGate(*swaps_arr[i]);
                *swapflags[i] = true;

              }

            }

            if(i % 2 == 0)
              contract_driver[driver_type]((void*)&tmpT1, Tout.first, (void*)&tmpT2, true);
            else
              contract_driver[driver_type](Tout.first, (void*)&tmpT1, (void*)&tmpT2, true);

          }

          if(i % 2 == 1){

            *(UniTensor<uni10_double64>*)Tout.first = tmpT1;

          }

        }else{

          driver_type = 3;

          UniTensor<uni10_complex128> tmpT1(*(UniTensor<uni10_complex128>*)tens_in_nd[0].first);
          UniTensor<uni10_complex128> tmpT2(*(UniTensor<uni10_complex128>*)tens_in_nd[1].first);

          if(tens_labels[0] == NULL){

            ((UniTensor<uni10_complex128>*)tens_in_nd[0].first)->clear();

          }
          else{

            tmpT1.setLabel(*tens_labels[0]);
            // AddSwaps
            uni10_error_msg(Qnum::isFermionic() && swapflags[0] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            uni10_error_msg(Qnum::isFermionic() && swaps_arr[0] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );

            if(Qnum::isFermionic() && !swapflags[0]){

              tmpT1.addGate(*swaps_arr[0]);
              *swapflags[0] = true;

            }

          }

          if(tens_labels[1] == NULL){

            ((UniTensor<uni10_complex128>*)tens_in_nd[1].first)->clear();

          }
          else{

            tmpT2.setLabel(*tens_labels[1]);
            // AddSwaps
            uni10_error_msg(Qnum::isFermionic() && swapflags[1] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            uni10_error_msg(Qnum::isFermionic() && swaps_arr[1] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
            if(Qnum::isFermionic() && !swapflags[1]){

              tmpT2.addGate(*swaps_arr[1]);
              *swapflags[1] = true;

            }

          }

          contract_driver[driver_type](Tout.first, (void*)&tmpT1, (void*)&tmpT2, true);

          uni10_uint64 i;

          for( i = 2; i < tens_in_nd.size(); i++){

            tmpT2 = *(UniTensor<uni10_complex128>*)tens_in_nd[i].first;

            if(tens_labels[i] == NULL){

              ((UniTensor<uni10_complex128>*)tens_in_nd[i].first)->clear();

            }
            else{

              tmpT2.setLabel(*tens_labels[i]);
              // AddSwaps
              uni10_error_msg(Qnum::isFermionic() && swapflags[i] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              uni10_error_msg(Qnum::isFermionic() && swaps_arr[i] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
              if(Qnum::isFermionic() && !swapflags[i]){

                tmpT2.addGate(*swaps_arr[i]);
                *swapflags[i] = true;

              }

            }

            if(i % 2 == 0)
              contract_driver[driver_type]((void*)&tmpT1, Tout.first, (void*)&tmpT2, true);
            else
              contract_driver[driver_type](Tout.first, (void*)&tmpT1, (void*)&tmpT2, true);

          }

          if(i % 2 == 1){

            *(UniTensor<uni10_complex128>*)Tout.first = tmpT1;

          }

        }

      }

    }

}; /* namespace uni10 */
