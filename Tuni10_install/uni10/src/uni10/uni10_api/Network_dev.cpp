/****************************************************************************
 *  @file Network_dev.cpp
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
#include <string>
#include <fstream>

#include "uni10/uni10_error.h"
#include "uni10/uni10_api/hirnk_linalg_inplace.h"
#include "uni10/uni10_api/Network_dev.h"


namespace uni10{

  std::ostream& operator<< (std::ostream& os, const Network_dev& net){

    if(!net.load){

      for(uni10_uint64 i = 0; i < net.names.size(); i++){
        os<<net.names[i]<< ": ";
        if(net.iBondNums[i])
          os<<"i[";
        for(int l = 0; l < net.iBondNums[i]; l++){
          os<<net.label_arr[i][l];
          if(l < net.iBondNums[i] - 1)
            os<<", ";

        } 
        if(net.iBondNums[i])
          os<<"] ";
        if(net.label_arr[i].size() - net.iBondNums[i])
          os<<"o[";
        for(uni10_uint64 l = net.iBondNums[i]; l < net.label_arr[i].size(); l++){
          os<<net.label_arr[i][l];
          if(l < net.label_arr[i].size() - 1)
            os<<", ";

        } 
        if(net.label_arr[i].size() - net.iBondNums[i])
          os<<"]";
        os<<std::endl;

      } 

    }
    else{

      std::vector<std::vector< std::string> > each_layer_out_strs(net.contraction_layers.size());
      std::map<std::string, std::map<uni10_int, uni10_uint64> > other_nd_labelDim;

      for(uni10_uint64 l = 0; l < net.contraction_layers.size(); l++)
        net.contraction_layers[l]->get_cout_strs(each_layer_out_strs[net.contraction_layers.size()-1-l], other_nd_labelDim);

      fprintf(stdout, "\n");

      for(uni10_uint64 ly = 0; ly < net.contraction_layers.size(); ly++){

        for(uni10_uint64 ln = 0; ln < each_layer_out_strs[ly].size(); ln++){

          fprintf(stdout, "%s\n", each_layer_out_strs[ly][ln].c_str());

        }

        fprintf(stdout, "\n");

      }

    }

    os << "";

    return os;

  }

  Network_dev::Network_dev(const std::string& _fname): fname(_fname), load(false), isChanged(false), ordered(false){

      fromfile(fname);

      uni10_int Tnum = label_arr.size() - 1;

      pos2pos_in_tree.assign(Tnum, std::vector<uni10_int>());

  }

  Network_dev::~Network_dev(){

    if(load)
      this->destruct();

  }

  void Network_dev::fromfile(const std::string& fname){  //names, name2pos, label_arr

    FILE *fp = fopen(fname.c_str(), "r");

    uni10_error_msg(!(fp), "Error in opening file ' %s '.", fname.c_str());

    int max_len = 1000;
    char tmp_line[max_len];
    char line[max_len];

    std::string name;
    std::vector<uni10_int> labels;

    int cnt = 0;

    while(fgets(tmp_line, max_len, fp)){

      uni10_error_msg((tmp_line[0] == '\t' || tmp_line[0] == '\n' || tmp_line[0] == ' '), "%s", "Can't have spaces, tabs or '\\n' at the head of line");

      if(tmp_line[0] == '#')
        continue;

      char *sharp_pch = strchr(tmp_line, '#');
      
      if(sharp_pch == NULL)
        strncpy(line, tmp_line, sizeof(tmp_line));
      else{
        strncpy(line, tmp_line, sharp_pch-tmp_line);
        line[sharp_pch-tmp_line] = '\0';
      }

      char* pch = strtok(line, ":");

      name.assign(pch);

      names.push_back(name);

      int iBonDNum = 0;
      bool in_iBD_rang = true;

      while(pch){

        pch = strtok(NULL, " , \n");

        if(pch){

          char* pch_end = pch+strlen(pch);
          char* semicolon_pos = std::find(pch, pch_end, ';');
          if( semicolon_pos != pch + strlen(pch)){

            int front_len = (semicolon_pos - pch);
            int back_len = (strlen(pch) - (semicolon_pos - pch) - 1);

            if(front_len > 0){
              labels.push_back(atoi(std::string(pch, semicolon_pos).c_str()));
              iBonDNum++;
            }

            if(back_len > 0){
              labels.push_back(atoi(std::string(semicolon_pos+1, pch_end).c_str()));
            }

            in_iBD_rang = false;

          }
          else
            labels.push_back(atoi(pch));

          if(in_iBD_rang)
            iBonDNum++;

        }

      }

      iBondNums.push_back(iBonDNum);

      name2pos[std::string(name)].push_back(cnt);

      label_arr.push_back(labels);

      labels.clear();
      // Check wether there is custimized order or not.
      if(name == "TOUT"){

        char order_line[1000];
        std::vector<std::string> order_dirty_str;

        if(fgets(order_line, 1000, fp)){

          //
          // remove (
          //
          for(uni10_uint64 c = 0; c < strlen(order_line); c++)
            if(order_line[c] == '(')
              order_line[c] = ' ';

          char *order_pch = strtok(order_line, " : \n");

          if(std::string(order_pch)== "ORDER"){

            while(order_pch){

              order_pch = strtok(NULL, "  \n");

              if(order_pch != NULL){
                order_dirty_str.push_back(order_pch);
              }

            }

            for(uni10_uint64 s = 0; s < order_dirty_str.size(); s++){

              uni10_uint64 sub_len = 0;

              std::string::iterator it = std::find(order_dirty_str[s].begin(), order_dirty_str[s].end(), ')');

              if(it == order_dirty_str[s].end()){

                uni10_error_msg(name2pos[order_dirty_str[s]].size() == 0, "There is no the tensor named \'%s\' in the netfile.", order_dirty_str[s].c_str());
                uni10_error_msg(name2pos[order_dirty_str[s]].size() != 1, "%s", "Please give the tensors different names if you want to set the contraction order manually.");
                contract_order.push_back(name2pos[order_dirty_str[s]][0]);

              }
              else{

                for(uni10_uint64 b = 0; b < order_dirty_str[s].size(); b++){

                  if(order_dirty_str[s][b] == ')'){

                    if(sub_len == 0){

                      contract_order.push_back(-1);

                    }else{

                      std::string Tname = order_dirty_str[s].substr(b-sub_len, sub_len);
                      uni10_error_msg(name2pos[Tname].size() != 1, "%s", "Please give the tensors different names if you want to set the contraction order manually.");
                      contract_order.push_back(name2pos[Tname][0]);

                      contract_order.push_back(-1);
                      sub_len = 0;

                    }

                  }
                  else
                    sub_len++;

                }

              }

            }

            ordered = true;
            break;

          }

        }

        break;

      }

      cnt++;

    } // End of loading net file.

    int numT = names.size() - 1;

    tensor_type.assign(numT, std::pair<void*, int>());

    std::vector<bool> found(numT, false);

    uni10_error_msg(!(names[numT] == "TOUT"), "Error in the network file ' %s '. Missing TOUT tensor. One must specify the bond labels for the resulting tensor by giving TOUT tag.\n  Hint: If the resulting tensor is a scalar(0-bond tensor), leave the TOUT empty like 'TOUT: '", fname.c_str());


    uni10_error_msg(!(names.size() > 2), "Error in the network file ' %s '. There must be at least two tensors in a tensor network.", fname.c_str());

    fclose(fp);
  }

  void Network_dev::init_layers(uni10_int& idx_ly, std::vector<uni10_int>& _contract_order, std::vector<uni10_int> &ori_ten_idx,
      std::vector<std::pair<const void*, int> >& _tens_type, std::vector<std::vector<uni10_int>* >& _label_arr, std::vector<std::string*>& _names){

    if(_tens_type.size() == 1)
      return ;

    std::vector<uni10_int> ly_nd_idx(3);

    std::vector<uni10_int> conOrder = _contract_order;

    //
    // Create a layer at the back of contration_layers.
    //
    contraction_layers.push_back((layer*)NULL);

    // Allocate for layer[idx_ly];
    contraction_layers[idx_ly] = new layer[1];

    layer* constructing_ly = contraction_layers[idx_ly];

    std::vector<uni10_int>                next_ly_contract_order;
    std::vector<std::pair<const void*, int> >   next_ly_tens_type;
    std::vector<std::vector<uni10_int>* > next_ly_label_arr;
    std::vector<std::string*>             next_ly_name;
    std::vector<uni10_int>                next_ly_ori_ten_idx;
    std::vector<uni10_int>                nd_pos_in_next_ly;

    // Push a dummy element back to avoid segmentation fault.
    conOrder.push_back(tensor_type.size());

    uni10_int num_of_tens_in_waiting_list = 0;
    uni10_int cnt_remaind_tens = 0;

    //
    uni10_int first_ten_idx_in_waiting_list = 0;
    uni10_int cnt_minus = 0;

    char tmp_nd_name[16];

    bool push_to_next_ly = false;

    for(uni10_int i = 0; i < (uni10_int)conOrder.size()-1; i++){

      if(conOrder[i] != -1){

        num_of_tens_in_waiting_list++;

      }else{

        cnt_minus++;

        //
        // There three conditions have to be considered if the next elemnet[i+1] is not -1, 
        // 1. num_of_ten_in_waiting_list - 1 >= cnt_minus;
        // 2. num_of_ten_in_waiting_list == 1 && cnt_minus == 1;
        // 3. num_of_ten_in_waiting_list - 1 < cnt_minus;
        //
        if(conOrder[i+1] != -1){

          // In the first case, 
          // 1. Create and push a node to the back of nds_in_layer in the constructing_ly.
          // 2. Assign N pair<void*, int>(), (vector<uni10_int>*)NULL and (std::sting*)NULL into the private vairables, tens_in_nd, *tens_labels and *tens_names, respectively.
          //    where N = cnt_minus + 1;
          // 3. Put (tensors, names, label_arr)[first_ten_idx_in_waiting_list + num_of_tens_in_waiting_list - 0],...
          //              ,[first_ten_idx_in_waiting_list + num_of_tens_in_waiting_list - cnt_minus],
          //    into the variables, tens_in_nd, tens_labels and tens_names, respectively.
          // 4. Consider the remaind idx, contract_order[first_ten_idx_in_wating_list + 0]...contract_order[first_ten_idx_in_wating_list + i- cnt_minus]. 
          //    We put the corresponded tensor, labels_arr and name into next_ly_tens_type, next_ly_label_arr and next_ly_name, respectively.
          // 5. Put Tout and name of this node to the back of next_ly_tens_type and next_ly_name, respectively.
          // 6. num_of_tens_in_waiting_list -> 0; 
          if(cnt_minus <= num_of_tens_in_waiting_list - 1 ){

            push_to_next_ly = false;

            constructing_ly->nds_in_layer.push_back((Node_dev*)NULL);
            constructing_ly->nds_in_layer.back() = new Node_dev[1];

            constructing_ly->nds_in_layer.back()->tens_in_nd.assign(cnt_minus+1, std::pair<const void*, int>());
            constructing_ly->nds_in_layer.back()->tens_labels.assign(cnt_minus+1, (std::vector<uni10_int>*)NULL);
            constructing_ly->nds_in_layer.back()->tens_names.assign(cnt_minus+1, (std::string*)NULL);

            uni10_int j = first_ten_idx_in_waiting_list;

            for(; j < i-(2*cnt_minus); j++){

              //std::cout << "J: " << j << std::endl;

              next_ly_contract_order.push_back(cnt_remaind_tens);
              cnt_remaind_tens++;

              next_ly_tens_type.push_back( _tens_type[ conOrder[j] ] );
              next_ly_label_arr.push_back( _label_arr[ conOrder[j] ] );
              next_ly_name.push_back( _names[ conOrder[j] ] );

              if(ori_ten_idx[ conOrder[j] ] != -1){

                next_ly_ori_ten_idx.push_back( ori_ten_idx[ conOrder[j] ] );

              }else{

                next_ly_ori_ten_idx.push_back( -1 );

              }

            }

            uni10_int cnt_tens_in_node = 0, end_idx = j+(cnt_minus) ;

            for(; j <  end_idx+1 ;j++){

              if(ori_ten_idx[ conOrder[end_idx - cnt_tens_in_node] ] != -1){

                ly_nd_idx[0] = idx_ly;
                ly_nd_idx[1] = constructing_ly->nds_in_layer.size() - 1;
                ly_nd_idx[2] = cnt_tens_in_node;
                this->pos2pos_in_tree[ ori_ten_idx[ conOrder[end_idx - cnt_tens_in_node] ] ] = ly_nd_idx;

              }

              constructing_ly->nds_in_layer.back()->tens_in_nd[cnt_tens_in_node] = _tens_type[ conOrder[end_idx - cnt_tens_in_node] ];
              constructing_ly->nds_in_layer.back()->tens_labels[cnt_tens_in_node] = _label_arr[ conOrder[end_idx - cnt_tens_in_node] ];
              constructing_ly->nds_in_layer.back()->tens_names[cnt_tens_in_node] = _names [ conOrder[end_idx - cnt_tens_in_node] ];
              cnt_tens_in_node++;

            }

            next_ly_contract_order.push_back(cnt_remaind_tens);
            nd_pos_in_next_ly.push_back(cnt_remaind_tens);
            cnt_remaind_tens++;

            next_ly_tens_type.push_back( constructing_ly->nds_in_layer.back()->Tout );
            next_ly_label_arr.push_back( (std::vector<uni10_int>*)NULL );

            sprintf(tmp_nd_name, "L%dN%ld", idx_ly, constructing_ly->nds_in_layer.size()-1);
            constructing_ly->nds_in_layer.back()->nd_name = std::string(tmp_nd_name);

            next_ly_name.push_back(&constructing_ly->nds_in_layer.back()->nd_name);
            next_ly_ori_ten_idx.push_back(-1);

            first_ten_idx_in_waiting_list = i+1;
            num_of_tens_in_waiting_list = 0;
            cnt_minus = 0;

          }

          else if(cnt_minus >= 1 && num_of_tens_in_waiting_list == 1 ){

            if(push_to_next_ly){

              next_ly_contract_order.push_back(cnt_remaind_tens);
              nd_pos_in_next_ly.push_back(cnt_remaind_tens);
              cnt_remaind_tens++;

              next_ly_tens_type.push_back( _tens_type[conOrder[first_ten_idx_in_waiting_list]]);
              next_ly_label_arr.push_back( _label_arr[conOrder[first_ten_idx_in_waiting_list]]);
              next_ly_name.push_back( _names[conOrder[first_ten_idx_in_waiting_list]]);

              uni10_int k;
              for(k = 0; k < cnt_minus; k++)
                next_ly_contract_order.push_back(-1);

              if(ori_ten_idx[ conOrder[first_ten_idx_in_waiting_list]] != -1){

                next_ly_ori_ten_idx.push_back( ori_ten_idx[ conOrder[first_ten_idx_in_waiting_list] ] );

              }else{

                next_ly_ori_ten_idx.push_back( -1 );

              }

            }else{

              if(ori_ten_idx[ conOrder[ first_ten_idx_in_waiting_list ] ] != -1){

                ly_nd_idx[0] = idx_ly;
                ly_nd_idx[1] = constructing_ly->nds_in_layer.size() - 1;
                ly_nd_idx[2] = constructing_ly->nds_in_layer.back()->tens_in_nd.size();
                this->pos2pos_in_tree[ ori_ten_idx[ conOrder[first_ten_idx_in_waiting_list] ] ] = ly_nd_idx;

              }

              constructing_ly->nds_in_layer.back()->tens_in_nd.push_back( _tens_type[ conOrder[first_ten_idx_in_waiting_list] ] );
              constructing_ly->nds_in_layer.back()->tens_labels.push_back( _label_arr[ conOrder[first_ten_idx_in_waiting_list] ] );
              constructing_ly->nds_in_layer.back()->tens_names.push_back( _names [ conOrder[first_ten_idx_in_waiting_list] ] );

              uni10_int k;
              for(k = 0; k < cnt_minus-num_of_tens_in_waiting_list; k++){
                next_ly_contract_order.push_back(-1);
              }

              if(cnt_minus > 1)
                push_to_next_ly = true;

            }

            first_ten_idx_in_waiting_list = i+1;
            num_of_tens_in_waiting_list = 0;
            cnt_minus = 0;

          }else if(cnt_minus > num_of_tens_in_waiting_list - 1 ){

            push_to_next_ly = false;

            constructing_ly->nds_in_layer.push_back((Node_dev*)NULL);
            constructing_ly->nds_in_layer.back() = new Node_dev[1];

            constructing_ly->nds_in_layer.back()->tens_in_nd.assign(num_of_tens_in_waiting_list, std::pair<const void*, int>());
            constructing_ly->nds_in_layer.back()->tens_labels.assign(num_of_tens_in_waiting_list, (std::vector<uni10_int>*)NULL);
            constructing_ly->nds_in_layer.back()->tens_names.assign(num_of_tens_in_waiting_list, (std::string*)NULL);

            uni10_int j = first_ten_idx_in_waiting_list;

            uni10_int cnt_tens_in_node = 0, end_idx = j + num_of_tens_in_waiting_list-1;

            for(; j < end_idx+1; j++){

              if(ori_ten_idx[ conOrder[ end_idx - cnt_tens_in_node ] ] != -1){

                ly_nd_idx[0] = idx_ly;
                ly_nd_idx[1] = constructing_ly->nds_in_layer.size() - 1;
                ly_nd_idx[2] = cnt_tens_in_node;
                this->pos2pos_in_tree[ ori_ten_idx[ conOrder[end_idx - cnt_tens_in_node] ] ] = ly_nd_idx;

              }

              constructing_ly->nds_in_layer.back()->tens_in_nd[cnt_tens_in_node] = _tens_type[ conOrder[end_idx - cnt_tens_in_node] ];
              constructing_ly->nds_in_layer.back()->tens_labels[cnt_tens_in_node] = _label_arr[ conOrder[end_idx - cnt_tens_in_node] ];
              constructing_ly->nds_in_layer.back()->tens_names[cnt_tens_in_node] = _names [ conOrder[end_idx - cnt_tens_in_node] ];
              cnt_tens_in_node++;

            }

            next_ly_contract_order.push_back(cnt_remaind_tens);
            nd_pos_in_next_ly.push_back(cnt_remaind_tens);
            cnt_remaind_tens++;

            next_ly_tens_type.push_back( constructing_ly->nds_in_layer.back()->Tout );
            next_ly_label_arr.push_back( (std::vector<uni10_int>*)NULL );

            sprintf(tmp_nd_name, "L%dN%ld", idx_ly, constructing_ly->nds_in_layer.size()-1);
            constructing_ly->nds_in_layer.back()->nd_name = std::string(tmp_nd_name);

            next_ly_name.push_back(&constructing_ly->nds_in_layer.back()->nd_name);
            next_ly_ori_ten_idx.push_back(-1);
     
            uni10_int k;

            for(k = 0; k < cnt_minus-num_of_tens_in_waiting_list+1; k++){

              next_ly_contract_order.push_back(-1);

            }

            first_ten_idx_in_waiting_list = i+1;
            num_of_tens_in_waiting_list = 0;
            cnt_minus = 0;

          }else{

            uni10_error_msg(true, "%s", "Unexpect errors. Please connect to the developers of uni10. Thanks");

          }

        }

      }

    }
 

    for(uni10_uint64 n = 0; n < constructing_ly->nds_in_layer.size(); n++){

      constructing_ly->nds_in_layer[n]->init_Tout();
      next_ly_tens_type[ nd_pos_in_next_ly[n] ] = constructing_ly->nds_in_layer[n]->Tout;

    }

    idx_ly++;

    init_layers(idx_ly, next_ly_contract_order, next_ly_ori_ten_idx, next_ly_tens_type, next_ly_label_arr, next_ly_name);


  }

  void Network_dev::pre_construct(bool force){

    if(load){

      if(isChanged){
        this->destruct_tree_only();
        isChanged=false;
        this->construct();
      }

      if(!force){
        this->destruct();
        this->construct();
      }

    }else{

      this->construct();

    }

  }

  void Network_dev::construct(){

    uni10_error_msg(load, "%s", "Unexpect error. Please connect to the developer of uni10 or create a issue on GitHub.");

    if(!ordered){

      NetOrder_dev opt_order(tensor_type, label_arr, names);
      contract_order = opt_order.generate_order();
      ordered = true;

    }
    
    uni10_int begin_ly = 0; 
    std::vector<std::vector<uni10_int>*> _label_arr;
    std::vector<std::string*> _names;
    std::vector<uni10_int> ori_ten_idx;

    for(uni10_uint64 i = 0; i < label_arr.size(); i++){
      ori_ten_idx.push_back(i);
      _label_arr.push_back(&label_arr[i]);
      _names.push_back(&names[i]);
    }

    // Init the contraction tree.
    init_layers(begin_ly, contract_order, ori_ten_idx, tensor_type, _label_arr, _names);
      
    if(Qnum::isFermionic())
      this->init_swap_arr_in_nodes();

    load = true;

  }

  void Network_dev::destruct(){

    //
    // Remove swspflags and swap_arr.
    //
    if(Qnum::isFermionic())
      this->del_swap_arr_in_nodes();
    //
    // Remove the contraction tree.
    //
    for(uni10_uint64 l = 0; l < contraction_layers.size(); l++){

      for(uni10_uint64 n = 0; n < contraction_layers[l]->nds_in_layer.size(); n++)
        delete [] contraction_layers[l]->nds_in_layer[n];

      delete [] contraction_layers[l];

    }

    contraction_layers.clear();

    load = false;

    //
    // Remove the contraction order.
    //
    contract_order.clear();
    ordered = false;

  }

  void Network_dev::destruct_tree_only(){

    //
    // Remove swspflags and swap_arr.
    //
    if(Qnum::isFermionic())
      this->del_swap_arr_in_nodes();
    //
    // Remove the contraction tree.
    //
    for(uni10_uint64 l = 0; l < contraction_layers.size(); l++){

      for(uni10_uint64 n = 0; n < contraction_layers[l]->nds_in_layer.size(); n++)
        delete [] contraction_layers[l]->nds_in_layer[n];

      delete [] contraction_layers[l];

    }

    contraction_layers.clear();

    load = false;

  }

  void Network_dev::putTensor_driver(int idx, const void* UniT, int typeID, bool force){

    uni10_error_msg(!(idx < (label_arr.size()-1)), "%s", "Index exceeds the number of the tensors in the list of network file.");

    if((!force) && load){
      destruct();
    }

    uni10_int RBondNum = (typeID == 1)? *(((UniTensor<uni10_double64>*)UniT)->RBondNum) : *(((UniTensor<uni10_complex128>*)UniT)->RBondNum);
    uni10_error_msg(!( RBondNum == iBondNums[idx] ), "The number of in-coming bonds does not match with the tensor ' %s ' specified in network file, %s.", names[idx].c_str(), fname.c_str());

    tensor_type[idx] = std::pair<const void*, int>(UniT, typeID);

    if(load){

      int ly_idx, nd_idx, ten_idx_in_nd;

      ly_idx        = pos2pos_in_tree[idx][0];
      nd_idx        = pos2pos_in_tree[idx][1];
      ten_idx_in_nd = pos2pos_in_tree[idx][2];

      uni10_int type_in_nd = contraction_layers[ly_idx]->nds_in_layer[nd_idx]->tens_in_nd[ten_idx_in_nd].second;

      if(tensor_type[idx].second != type_in_nd)
        isChanged = true;

      contraction_layers[ly_idx]->nds_in_layer[nd_idx]->tens_in_nd[ten_idx_in_nd] = tensor_type[idx];

      if(Qnum::isFermionic()){

        bool *fptr = contraction_layers[ly_idx]->nds_in_layer[nd_idx]->swapflags[ten_idx_in_nd];
        std::vector<_Swap> *arrptr= contraction_layers[ly_idx]->nds_in_layer[nd_idx]->swaps_arr[ten_idx_in_nd];

        uni10_error_msg(fptr==NULL || arrptr==NULL, "%s", "Unexpected error. Please connect to the developers of uni10.");

        *fptr = false;

      }

    }

  }

  void Network_dev::putTensor(int idx, const UniTensor<uni10_double64>& UniT, bool force){

    this->putTensor_driver(idx, (void*)&UniT, UniT.typeID(), force);

  }

  void Network_dev::putTensor(int idx, const UniTensor<uni10_complex128>& UniT, bool force){

    this->putTensor_driver(idx, (void*)&UniT, UniT.typeID(), force);

  }

  void Network_dev::putTensor(const std::string& name, const UniTensor<uni10_double64>& UniT, bool force){

    std::map<std::string, std::vector<uni10_int> >::const_iterator it = name2pos.find(name);
    uni10_error_msg(!(it != name2pos.end()), "There is no tensor named ' %s ' in the network file", name.c_str());

    uni10_uint64 i = 0;
    for(; i < it->second.size(); i++)
      putTensor(it->second[i], UniT, force);

  }

  void Network_dev::putTensor(const std::string& name, const UniTensor<uni10_complex128>& UniT, bool force){

    std::map<std::string, std::vector<uni10_int> >::const_iterator it = name2pos.find(name);
    uni10_error_msg(!(it != name2pos.end()), "There is no tensor named ' %s ' in the network file", name.c_str());

    uni10_uint64 i = 0;

    for(; i < it->second.size(); i++)
      putTensor(it->second[i], UniT, force);

  }

  void Network_dev::putTensorT(int idx, const UniTensor<uni10_double64>& UniT, bool force){

    UniTensor<uni10_double64> transT;
    transpose(transT, UniT, INPLACE) ;
    putTensor(idx, transT, force);

  }

  void Network_dev::putTensorT(int idx, const UniTensor<uni10_complex128>& UniT, bool force){

    UniTensor<uni10_complex128> transT;
    transpose(transT, UniT, INPLACE) ;
    putTensor(idx, transT, force);

  }

  void Network_dev::putTensorT(const std::string& nameT, const UniTensor<uni10_double64>& UniT, bool force){

    UniTensor<uni10_double64> transT;
    transpose(transT, UniT, INPLACE) ;
    putTensor(nameT, transT, force);

  }

  void Network_dev::putTensorT(const std::string& nameT, const UniTensor<uni10_complex128>& UniT, bool force){

    UniTensor<uni10_complex128> transT;
    transpose(transT, UniT, INPLACE) ;
    putTensor(nameT, transT, force);

  }

  void Network_dev::putTensorD(int idx, const UniTensor<uni10_double64>& UniT, bool force){

    UniTensor<uni10_double64> transT;
    dagger(transT, UniT, INPLACE) ;
    putTensor(idx, transT, force);


  }

  void Network_dev::putTensorD(int idx, const UniTensor<uni10_complex128>& UniT, bool force){

    UniTensor<uni10_complex128> transT;
    dagger(transT, UniT, INPLACE) ;
    putTensor(idx, transT, force);

  }

  void Network_dev::putTensorD(const std::string& nameT, const UniTensor<uni10_double64>& UniT, bool force){

    UniTensor<uni10_double64> transT;
    dagger(transT, UniT, INPLACE) ;
    putTensor(nameT, transT, force);

  }

  void Network_dev::putTensorD(const std::string& nameT, const UniTensor<uni10_complex128>& UniT, bool force){

    UniTensor<uni10_complex128> transT;
    dagger(transT, UniT, INPLACE) ;
    putTensor(nameT, transT, force);

  }

  void Network_dev::launch(UniTensor<uni10_double64>& Tout, const std::string& Tname){

    if(isChanged && load){
      this->destruct_tree_only();
      isChanged = false;
    }

    if(!load) 
      this->construct();

    uni10_uint64 i;

    for(i = 0; i < contraction_layers.size(); i++){
      contraction_layers[i]->merge();
    }

    UniTensor<uni10_double64> tmp(*(UniTensor<uni10_double64>*)contraction_layers.back()->nds_in_layer[0]->Tout.first);
    permute(Tout, tmp, label_arr.back(), iBondNums.back(), INPLACE);
    Tout.setName(Tname);

    ((UniTensor<uni10_double64>*)contraction_layers.back()->nds_in_layer[0]->Tout.first)->clear();

  }

  void Network_dev::launch(UniTensor<uni10_complex128>& Tout, const std::string& Tname){

    if(isChanged && load){
      this->destruct_tree_only();
      isChanged = false;
    }

    if(!load) 
      this->construct();

    uni10_uint64 i;

    for(i = 0; i < contraction_layers.size(); i++)
      contraction_layers[i]->merge();

    UniTensor<uni10_complex128> tmp(*(UniTensor<uni10_complex128>*)contraction_layers.back()->nds_in_layer[0]->Tout.first);
    permute(Tout, tmp, label_arr.back(), iBondNums.back(), INPLACE);
    Tout.setName(Tname);

    ((UniTensor<uni10_complex128>*)contraction_layers.back()->nds_in_layer[0]->Tout.first)->clear();

  }

  std::string Network_dev::get_contract_order() const{

    std::string contract_order_str = names[contract_order[0]] + " ";

    for(uni10_uint64 i = 0; i < contract_order.size()-1; i++){

      if(contract_order[i+1] == -1)
        contract_order_str += "# ";
      else
        contract_order_str += (names[contract_order[i+1]] + " ");
    }

    return contract_order_str;

  }

  //
  // Before initialize the swap_arr in each nodes, 
  // the contraction_layers should be build and obtain the value of pos2pos_in_tree firstly.
  // 
  void Network_dev::init_swap_arr_in_nodes(){

    for(uni10_uint64 l = 0; l < contraction_layers.size(); l++){

      layer *lptr = contraction_layers[l];

      for(uni10_uint64 n = 0; n < lptr->nds_in_layer.size(); l++){

        Node_dev *ndptr = lptr->nds_in_layer[n];

        uni10_error_msg(ndptr->swapflags.size() != 0, "%s", "Unexpected error. Please connect to the developers of uni10." );
        uni10_error_msg(ndptr->swaps_arr.size() != 0, "%s", "Unexpected error. Please connect to the developers of uni10." );
        
        ndptr->swapflags.assign(ndptr->tens_in_nd.size(), (bool*)NULL);
        ndptr->swaps_arr.assign(ndptr->tens_in_nd.size(), (std::vector<_Swap>*)NULL);

      }

    }

    for(uni10_uint64 i = 0; i < pos2pos_in_tree.size(); i++){

      layer    *lptr  = contraction_layers[pos2pos_in_tree[i][0]];
      Node_dev *ndptr = lptr->nds_in_layer[pos2pos_in_tree[i][1]];
  
      uni10_error_msg(ndptr->tens_labels[pos2pos_in_tree[i][2]] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );

      ndptr->swapflags[pos2pos_in_tree[i][2]] = new bool[1];
     *ndptr->swapflags[pos2pos_in_tree[i][2]] = false;
      ndptr->swaps_arr[pos2pos_in_tree[i][2]] = new std::vector<_Swap>[1];

      // Then use addSwap to assign the _Swaps into swaps_arr.
      this->addSwap();

    }

  }

  void Network_dev::del_swap_arr_in_nodes(){

    uni10_error_msg(pos2pos_in_tree.size() != label_arr.size()-1, "%s", "Unexpected error. Please connect to the developers of uni10.");

    for(uni10_uint64 i = 0; i < pos2pos_in_tree.size(); i++){

      layer    *lptr  = contraction_layers[pos2pos_in_tree[i][0]];
      Node_dev *ndptr = lptr->nds_in_layer[pos2pos_in_tree[i][1]];

      uni10_error_msg(ndptr->swapflags[pos2pos_in_tree[i][2]] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );
      uni10_error_msg(ndptr->swaps_arr[pos2pos_in_tree[i][2]] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );

      delete [] ndptr->swapflags[pos2pos_in_tree[i][2]];
      delete [] ndptr->swaps_arr[pos2pos_in_tree[i][2]];

    }

    for(uni10_uint64 l = 0; l < contraction_layers.size(); l++){

      layer *lptr = contraction_layers[l];

      for(uni10_uint64 n = 0; n < lptr->nds_in_layer.size(); l++){

        Node_dev *ndptr = lptr->nds_in_layer[n];

        ndptr->swapflags.clear();
        ndptr->swaps_arr.clear();

      }

    }

  }

  //
  // Before using this function to add _Swap into a specified noda, 
  // we should init_layers and init_swap_arr in the node at first.
  //
  void Network_dev::addSwap(){

    uni10_error_msg(!load, "%s", "Unexpected error. Please connect to the developers of uni10.");

    int Tnum = label_arr.size()-1;

    std::vector< std::vector<_Swap> > swaps_arr(Tnum);

    std::vector<uni10_int> tenOrder(Tnum);

    uni10_uint64 cnt = 0;

    for(uni10_uint64 i = 0; i < contract_order.size(); i++){

      if(contract_order[i] != -1){
        tenOrder[cnt] = contract_order[i];
        cnt++;
      }

    }

    std::vector<_Swap> tenSwaps = recSwap(tenOrder);

    std::vector<_Swap> swtmp;

    for(int s = 0; s < (int)tenSwaps.size(); s++){

      int driver_type = (tensor_type[tenSwaps[s].b1].second - 1) + (2 * (tensor_type[tenSwaps[s].b2].second - 1) );

      if(driver_type == 0){

        swtmp = ((UniTensor<uni10_double64>*)tensor_type[tenSwaps[s].b1].first)->exSwap(*((UniTensor<uni10_double64>*)tensor_type[tenSwaps[s].b2].first));

      }
      else if(driver_type == 1){

        swtmp = ((UniTensor<uni10_complex128>*)tensor_type[tenSwaps[s].b1].first)->exSwap(*((UniTensor<uni10_double64>*)tensor_type[tenSwaps[s].b2].first));

      }
      else if(driver_type == 2){

        swtmp = ((UniTensor<uni10_double64>*)tensor_type[tenSwaps[s].b1].first)->exSwap(*((UniTensor<uni10_complex128>*)tensor_type[tenSwaps[s].b2].first));
      }

      else if(driver_type == 3){

        swtmp = ((UniTensor<uni10_complex128>*)tensor_type[tenSwaps[s].b1].first)->exSwap(*((UniTensor<uni10_complex128>*)tensor_type[tenSwaps[s].b2].first));
      }
      else{

        uni10_error_msg(true, "%s", "Unexpected error. Please, connect to the developers of uni10.")

      }

      swaps_arr[tenSwaps[s].b1].insert(swaps_arr[tenSwaps[s].b1].end(), swtmp.begin(), swtmp.end());

    }

    //Distinct Swaps of each tensors
    
    for(uni10_int t = 0; t < Tnum; t++){

      std::map<uni10_int, bool> recs;
      std::map<uni10_int, bool>::iterator it;

      
      uni10_uint64 bondNum = tensor_type[t].second == 1 ? 
        ((UniTensor<uni10_double64>*)tensor_type[t].first)->bonds->size() : ((UniTensor<uni10_complex128>*)tensor_type[t].first)->bonds->size();

      uni10_int is, ib;
      //int hash;
      for(uni10_int s = 0; s < (uni10_int)swaps_arr[t].size(); s++){
        if(swaps_arr[t][s].b1 < swaps_arr[t][s].b2){
          is = swaps_arr[t][s].b1;
          ib = swaps_arr[t][s].b2;
        }
        else{
          is = swaps_arr[t][s].b2;
          ib = swaps_arr[t][s].b1;
        }
        uni10_int hash = is * bondNum + ib;
        if((it = recs.find(hash)) != recs.end())
          it->second ^= true;
        else
          recs[hash] = true;
      }

      swaps_arr[t].clear();
      _Swap sp;

      for (it=recs.begin(); it!=recs.end(); it++){
        if(it->second){
          sp.b1 = (it->first) / bondNum;
          sp.b2 = (it->first) % bondNum;
          swaps_arr[t].push_back(sp);
        }
      }
    }

    for(uni10_int i = 0; i < (uni10_int)swaps_arr.size(); i++){
     
      layer    *lptr  = contraction_layers[pos2pos_in_tree[i][0]];
      Node_dev *ndptr = lptr->nds_in_layer[pos2pos_in_tree[i][1]];
  
      uni10_error_msg(ndptr->swaps_arr[pos2pos_in_tree[i][2]] == NULL, "%s", "Unexpected error. Please connect to the developers of uni10." );

      *ndptr->swaps_arr[pos2pos_in_tree[i][2]] = swaps_arr[i];

    }

  }


}; /* namespace uni10 */
