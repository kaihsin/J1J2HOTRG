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
 *  @author Yun-Da Hsieh, Ying-Jer Kao
 *  @date 2014-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/
#include <string>
#include <fstream>

#include "uni10/uni10_error.h"
#include "uni10/uni10_api/hirnk_linalg.h"
#include "uni10/uni10_api/Network.h"
#include "uni10/uni10_api/network_tools/network_tools.h"
#include "uni10/uni10_api/network_tools/NetOrder.h"


namespace uni10{

  template <typename uni10_type>
    Network<uni10_type>::Network(const std::string& _fname): fname(_fname), root(NULL), load(false), times(0), tot_elem(0), max_elem(0), hasOrder(false), order_pos(0){
      fromfile(fname);
      int Tnum = label_arr.size() - 1;
      swapflags.assign(Tnum, false);
      std::vector<_Swap> swaps;
      swaps_arr.assign(Tnum, swaps);
      leafs.assign(Tnum, (Node<uni10_type>*)NULL);
      tensors.assign(Tnum, (UniTensor<uni10_type>*)NULL);
    }

  template <typename uni10_type>
    Network<uni10_type>::Network(const std::string& _fname, const std::vector<UniTensor<uni10_type>*>& tens): fname(_fname), root(NULL), load(false), times(0), tot_elem(0), max_elem(0), hasOrder(false), order_pos(0){
      fromfile(fname);
      uni10_error_msg(!((label_arr.size() - 1) == tens.size()), "The size of the input vector of tensors does not match for the number of tensors in the network file ' %s '.", fname.c_str());
      int Tnum = tens.size();
      swapflags.assign(Tnum, false);
      std::vector<_Swap> swaps;
      swaps_arr.assign(Tnum, swaps);
      leafs.assign(Tnum, (Node<uni10_type>*)NULL);
      tensors.assign(Tnum, (UniTensor<uni10_type>*)NULL);
      for(int i = 0; i < Tnum; i++){
        uni10_error_msg(!(*tens[i]->RBondNum == Rnums[i]), "The number of in-coming bonds does not match with the tensor ' %s ' specified in network file ' %s '.", names[i].c_str(), fname.c_str());
        UniTensor<uni10_type>* ten = new UniTensor<uni10_type>(*(tens[i]));
        ten->setName(names[i]);
        ten->setLabel(label_arr[i]);
        tensors[i] = ten;
        Node<uni10_type>* ndp = new Node<uni10_type>(ten);
        leafs[i] = ndp;
      }

      if(!hasOrder){

        names = std::vector<std::string>();
        label_arr = std::vector< std::vector<int> >();

        FILE* fp = fopen(fname.c_str(), "rb+");
        NetOrder _order(tensors);
        char* netorder = _order.generate_order();
        fseek(fp, order_pos, SEEK_SET);
        fprintf(fp, "ORDER: %s", netorder);
        fclose(fp);

        fromfile(fname);

      }

      construct();
    }

  template <typename uni10_type>
    void Network<uni10_type>::fromfile(const std::string& fname){//names, name2pos, label_arr, Rnums, order, brakets
      std::string str;
      std::ifstream infile;
      infile.open (fname.c_str());
      uni10_error_msg(!(infile.is_open()), "Error in opening file ' %s '.", fname.c_str());
      int lnum = 0;
      int MAXLINES = 1000;
      int pos = 0;
      int endpos = 0;
      std::string tar("1234567890-");
      std::vector<std::string> ord;
      while(lnum < MAXLINES){
        getline(infile, str); // Saves the line in STRING.
        endpos = 0;
        if(infile.eof())
          break;
        pos = str.find(":");
        if(pos == (int)std::string::npos)
          break;
        std::string name = str.substr(0, pos);
      trim(name);
      if(name == "ORDER"){
        std::string bra("(");
        if(str.find(bra, pos+1) == std::string::npos){
          std::string del(" ,;");
          while(((pos = str.find_first_not_of(del, pos + 1)) != (int)std::string::npos)){
            endpos = str.find_first_of(del, pos + 1);
            if(endpos == (int)std::string::npos)
              ord.push_back(str.substr(pos));
            else
              ord.push_back(str.substr(pos, endpos - pos));
            pos = endpos;
            if(pos == (int)std::string::npos)
              break;
          }
        }
        else{
          std::string del(" ,;*()");
          endpos = str.find_first_of(del, pos + 1);
          while(((pos = str.find_first_not_of(del, pos + 1)) != (int)std::string::npos)){
            std::string bras = str.substr(endpos, pos - endpos);
            int minus = std::count(bras.begin(), bras.end(), ')');
            int plus = std::count(bras.begin(), bras.end(), '(');
            if(minus)
              brakets.push_back(-minus);
            if(plus)
              brakets.push_back(plus);
            endpos = str.find_first_of(del, pos + 1);
            if(endpos == (int)std::string::npos){
              ord.push_back(str.substr(pos));
              brakets.push_back(0);
            }
            else{
              ord.push_back(str.substr(pos, endpos - pos));
              brakets.push_back(0);
            }
            pos = endpos;
            if(pos == (int)std::string::npos)
              break;
          }
          if(endpos != (int)std::string::npos){
            std::string bras = str.substr(endpos);
            int minus = std::count(bras.begin(), bras.end(), ')');
            if(minus)
              brakets.push_back(-minus);
          }
          int sum = 0;
          for(int i = 0; i < (int)brakets.size(); i++){
            sum += brakets[i];
          }

          uni10_error_msg(!(sum == 0), "Error in the network file ' %s '. There are imbalance brackets when specifying the contraction order.", fname.c_str());
        }
        hasOrder = true;
        break;
      }
      name2pos[name] = names.size();
      if(name == "TOUT")
        order_pos = infile.tellg();
      names.push_back(name);
      std::vector<int> labels;
      int Rnum = -1;
      int cnt = 0;
      int tmp;
      while((pos = str.find_first_of(tar, pos + 1)) != (int)std::string::npos){
        if(Rnum == -1){
          tmp = str.find(";", endpos);
          if(tmp != (int)std::string::npos && tmp < pos)
            Rnum = cnt;
        }
        endpos = str.find_first_not_of(tar, pos + 1);
        std::string label;
        if(endpos == (int)std::string::npos)
          label = str.substr(pos);
        else
          label = str.substr(pos, endpos - pos);
        char* pEnd;
        labels.push_back(strtol(label.c_str(), &pEnd, 10));
        pos = endpos;
        if(Rnum == -1)
          cnt++;
        if(pos == (int)std::string::npos)
          break;
      }
      label_arr.push_back(labels);
      if(Rnum == -1)
        Rnum = labels.size();
      Rnums.push_back(Rnum);
      lnum ++;
      }
      int numT = names.size() - 1;
      std::vector<bool> found(numT, false);

      uni10_error_msg(!(names[numT] == "TOUT"), "Error in the network file ' %s '. Missing TOUT tensor. One must specify the bond labels for the resulting tensor by giving TOUT tag.\n  Hint: If the resulting tensor is a scalar(0-bond tensor), leave the TOUT empty like 'TOUT: '", fname.c_str());

      uni10_error_msg(!(names.size() > 2), "Error in the network file ' %s '. There must be at least two tensors in a tensor network.", fname.c_str());
      order.assign(numT, 0);
      if(ord.size() > 0){
        uni10_error_msg(!((int)ord.size() == numT), "Error in the network file ' %s '. Some tensors are missing in the contraction order.", fname.c_str());
      std::map<std::string, size_t>::iterator it;
      for(int i = 0; i < numT; i++){
        it = name2pos.find(ord[i]);
        uni10_error_msg(!(it != name2pos.end()), "Error in the network file ' %s '. ' %s ' in the contraction order is not in the list of tensors above.", fname.c_str(), ord[i].c_str());
        order[i] = it->second;
        uni10_error_msg(!(found[order[i]] == false), "Error in the network file ' %s '. ' %s ' appears more than once in the contraction order.", fname.c_str(), ord[i].c_str());
        found[order[i]] = true;
      }
      }
      else{
        for(int i = 0; i < numT; i++)
          order[i] = i;
      }
      infile.close();
    }

  template <typename uni10_type>
    void Network<uni10_type>::construct(){
      if(brakets.size()){
        std::vector<Node<uni10_type>*> stack(leafs.size(), NULL);
        int cursor = 0;
        int cnt = 0;
        for(int i = 0; i < (int)brakets.size(); i++){
          if(brakets[i] < 0){
            for(int c = 0; c < -brakets[i]; c++){
              Node<uni10_type>* par = new Node<uni10_type>(stack[cursor - 2]->contract(stack[cursor - 1]));
              par->left = stack[cursor - 2];
              par->right = stack[cursor - 1];
              par->left->parent = par;
              par->right->parent = par;
              stack[cursor - 2] = par;
              cursor--;
            if(cursor < 2)	//prevent breakdown because of extra braket
              break;
          }
        }
        else if(brakets[i] == 0){
          uni10_error_msg(leafs[order[cnt]] == NULL, "(((Tensor ' %s 'has not yet been given.\n  Hint: Use addTensor() to add a tensor to a network.\n", names[order[cnt]].c_str());
          stack[cursor] = leafs[order[cnt]];
          cnt++;
          cursor++;
        }
      }
      while(cursor > 1){//for imcomplete brakets
        Node<uni10_type>* par = new Node<uni10_type>(stack[cursor - 2]->contract(stack[cursor - 1]));
        par->left = stack[cursor - 2];
        par->right = stack[cursor - 1];
        par->left->parent = par;
        par->right->parent = par;
        stack[cursor - 2] = par;
        cursor--;
      }
      root = stack[0];
    }
    else{
      for(int i = 0; i < (int)order.size(); i++){
        uni10_error_msg(leafs[order[i]] == NULL, "Tensor ' %s ' has not yet been given.\n  Hint: Use putTensor() to add a tensor to a network.\n", names[order[i]].c_str());
        matching(leafs[order[i]], root);
      }
    }
    int Tnum = label_arr.size() - 1;
    if(root->labels.size() == label_arr[Tnum].size()){
      for(int l = 0; l < (int)root->labels.size(); l++){
        bool found = false;
        for(int t = 0; t < (int)label_arr[Tnum].size(); t++)
          if(root->labels[l] == label_arr[Tnum][t]){
            found = true;
            break;
          }
        if(!found){
          char err[2048];
          char buf[8];
          sprintf(err, "Error when constructing the network. The labels of the resulting tensor, ( ");
          for(int i = 0; i < (int)root->labels.size(); i++){
            sprintf(buf, "%d ", root->labels[i]);
            strcat(err, buf);
          }
          strcat(err,"), do not match with the labels of 'TOUT' in the network file");
          uni10_error_msg(true, "%s", err);
        }
      }

    }
    else{
      uni10_error_msg(true, "%s", "Error when constructing the network. The bond number of the resulting tensor is different from the bond number of 'TOUT'");
    }
    addSwap();
    load = true;

    }

  template <typename uni10_type>
    void Network<uni10_type>::putTensor(size_t idx, const UniTensor<uni10_type>* UniT, bool force){

      uni10_error_msg(!(idx < (label_arr.size()-1)), "%s", "Index exceeds the number of the tensors in the list of network file.");

      if((!force) && load){
        destruct();
      }

      uni10_error_msg(!(*UniT->RBondNum == Rnums[idx]), "The number of in-coming bonds does not match with the tensor ' %s ' specified in network file", names[idx].c_str());

      if(leafs[idx] != NULL){
        *(tensors[idx]) = *UniT;
        tensors[idx]->setLabel(label_arr[idx]);
        tensors[idx]->setName(names[idx]);
        swapflags[idx] = false;
      }
      else{
        UniTensor<uni10_type>* ten = new UniTensor<uni10_type>(*UniT);
        ten->setName(names[idx]);
        ten->setLabel(label_arr[idx]);
        tensors[idx] = ten;
        Node<uni10_type>* ndp = new Node<uni10_type>(ten);
        leafs[idx] = ndp;
      }
    }

  template <typename uni10_type>
    void Network<uni10_type>::putTensor(size_t idx, const UniTensor<uni10_type>& UniT, bool force){
      putTensor(idx, &UniT, force);
    }

  template <typename uni10_type>
    void Network<uni10_type>::putTensor(const std::string& name, const UniTensor<uni10_type>* UniT, bool force){

      std::map<std::string, size_t>::const_iterator it = name2pos.find(name);
      uni10_error_msg(!(it != name2pos.end()), "There is no tensor named ' %s ' in the network file", name.c_str());
      putTensor(it->second, UniT, force);

    }

  template <typename uni10_type>
    void Network<uni10_type>::putTensor(const std::string& name, const UniTensor<uni10_type>& UniT, bool force){
      putTensor(name, &UniT, force);
    }

  template <typename uni10_type>
    void Network<uni10_type>::putTensorT(const std::string& nameT, const UniTensor<uni10_type>* UniT, bool force){
      std::map<std::string, size_t>::const_iterator itT = name2pos.find(nameT);
      uni10_error_msg(!(itT != name2pos.end()), "There is no tensor named ' %s ' in the network file", nameT.c_str());
      UniTensor<uni10_type> transT;
      transpose(transT, *UniT, INPLACE) ;
      putTensor(itT->second, &transT, force);
    }

  template <typename uni10_type>
    void Network<uni10_type>::putTensorT(const std::string& nameT, const UniTensor<uni10_type>& UniT, bool force){
      putTensorT(nameT, &UniT, force);
    }

  template <typename uni10_type>
    void Network<uni10_type>::branch(Node<uni10_type>* sbj, Node<uni10_type>* tar){
      Node<uni10_type>* par = new Node<uni10_type>(tar->contract(sbj));
      if(sbj->parent == NULL){  //create a parent node
        if(tar->parent != NULL){  //tar is not root
          if(tar->parent->left == tar)  // tar on the left side of its parent
            tar->parent->left = par;
          else
            tar->parent->right = par;
          par->parent = tar->parent;
        }
        else{ //tar is root
          par->parent = NULL;
          root = par;
        }
      }
      else{ //sbj and tar have same parent and replace the parent node
        if(tar->parent->parent != NULL){
          if(tar->parent->parent->left == tar->parent)  // tar on the left side of its parent
            tar->parent->parent->left = par;
          else
            tar->parent->parent->right = par;
          par->parent = tar->parent->parent;
        }
        else{ //tar->parent is root
          par->parent = NULL;
          root = par;
        }
        delete tar->parent;
      }
      par->left = tar;
      par->right = sbj;
      tar->parent = par;
      sbj->parent = par;
      par->point = tar->metric(sbj);
      if(sbj->parent->parent != NULL){  //propagate up
        sbj = sbj->parent;
        branch(sbj->parent->right, sbj->parent->left);
      }
    }

  template <typename uni10_type>
    void Network<uni10_type>::matching(Node<uni10_type>* sbj, Node<uni10_type>* tar){
      if(tar == NULL){  //tar is root
        root = sbj;
      }
      else if(tar->T == NULL){  //not leaf
        if(sbj->metric(tar) > 0){ //has contracted bonds

          uni10_error_msg(!(tar->left != NULL && tar->right != NULL), "%s", "Fatal error(code = N1). Please contact the developer of the uni10 library.");

          uni10_float32 tar_p = tar->point;
          uni10_float32 lft_p = 0, rht_p = 0;
          if((lft_p = sbj->metric(tar->left)) > tar_p || (rht_p = sbj->metric(tar->right)) > tar_p){	//go deeper comparison to the children
            if(lft_p > rht_p)
              matching(sbj, tar->left);
            else
              matching(sbj, tar->right);
          }
          else  //contract
            branch(sbj, tar);
        }
        else  //contract
          branch(sbj, tar);
      }
      else{ //contract
        branch(sbj, tar);
      }
    }

  template <typename uni10_type>
    void Network<uni10_type>::clean(Node<uni10_type>* nd){
      if(nd->T != NULL)	//leaf
        return;
      clean(nd->left);
      clean(nd->right);
      delete nd;
    }

  template <typename uni10_type>
    void Network<uni10_type>::destruct(){
      clean(root);
      root = NULL;
      for(int i = 0; i < (int)leafs.size(); i++)
        leafs[i]->delink();
      conOrder.clear();
      for(int t = 0; t < (int)tensors.size(); t++){
        if(Qnum::isFermionic() && swapflags[t]){
          tensors[t]->addGate(swaps_arr[t]);
          swapflags[t] = false;
        }
        swaps_arr[t].clear();
      }
      load = false;
    }


  template <typename uni10_type>
    UniTensor<uni10_type> Network<uni10_type>::launch(const std::string& _name){

      if(!hasOrder){

        names = std::vector<std::string>();
        label_arr = std::vector< std::vector<int> >();

        FILE* fp = fopen(fname.c_str(), "rb+");
        NetOrder _order(tensors);
        char* netorder = _order.generate_order();
        fseek(fp, order_pos, SEEK_SET);
        fprintf(fp, "ORDER: %s\n", netorder);
        fclose(fp);

        fromfile(fname);

      }

      if(!load)
        construct();

      for(int t = 0; t < (int)tensors.size(); t++)
        if(Qnum::isFermionic() && !swapflags[t]){
          tensors[t]->addGate(swaps_arr[t]);
          swapflags[t] = true;
        }
      UniTensor<uni10_type> UniT = merge(root);
      int idx = label_arr.size() - 1;
      if(label_arr.size() > 0 && label_arr[idx].size() > 1)
        UniT = permute(UniT, label_arr[idx], Rnums[idx]);
      UniT.setName(_name);
      return UniT;


    }

  template <typename uni10_type>
   void Network<uni10_type>::launch(UniTensor<uni10_type>& Tout, const std::string& _name){

      if(!hasOrder){

        names = std::vector<std::string>();
        label_arr = std::vector< std::vector<int> >();

        FILE* fp = fopen(fname.c_str(), "rb+");
        NetOrder _order(tensors);
        char* netorder = _order.generate_order();
        fseek(fp, order_pos, SEEK_SET);
        fprintf(fp, "ORDER: %s\n", netorder);
        fclose(fp);

        fromfile(fname);

      }

      if(!load)
        construct();

      for(int t = 0; t < (int)tensors.size(); t++)
        if(Qnum::isFermionic() && !swapflags[t]){
          tensors[t]->addGate(swaps_arr[t]);
          swapflags[t] = true;
        }
      UniTensor<uni10_type> UniT = merge(root);
      int idx = label_arr.size() - 1;
      if(label_arr.size() > 0 && label_arr[idx].size() > 1)
        permute(Tout, UniT, label_arr[idx], Rnums[idx], INPLACE);
      else
        Tout = UniT;
      Tout.setName(_name);

    }

  template <typename uni10_type>
    UniTensor<uni10_type> Network<uni10_type>::merge(Node<uni10_type>* nd){
      if(nd->left->T == NULL){
        UniTensor<uni10_type> lftT = merge(nd->left);
        if(nd->right->T == NULL){
          UniTensor<uni10_type> rhtT = merge(nd->right);
          return contract(lftT, rhtT, true);
        }
        else{
          return contract(lftT, *(nd->right->T), true);
        }
      }
      else{
        if(nd->right->T == NULL){
          UniTensor<uni10_type> rhtT = merge(nd->right);
          return contract(*(nd->left->T), rhtT, true);
        }
        else{
          return contract(*(nd->left->T), *(nd->right->T), true);
        }
      }
    }

  template <typename uni10_type>
    Network<uni10_type>::~Network(){
      if(load)
        destruct();
      for(int i = 0; i < (int)leafs.size(); i++)
        delete leafs[i];
      for(int i = 0; i < (int)tensors.size(); i++)
        delete tensors[i];
    }

  template <typename uni10_type>
    int Network<uni10_type>::rollcall(){
      if(!load){
        for(int i = 0; i < (int)leafs.size(); i++)
          if(leafs[i] == NULL){
            return i;
          }
        construct();
      }
      return -1;
    }

  template <typename uni10_type>
    size_t Network<uni10_type>::sum_of_memory_usage(){
      if(rollcall() >= 0)
        return 0;
      return _sum_of_tensor_elem(root) * sizeof(uni10_type);
    }

  template <typename uni10_type>
    size_t Network<uni10_type>::_sum_of_tensor_elem(Node<uni10_type>* nd) const{
      if(nd == NULL)
        return 0;
      return nd->elemNum + _sum_of_tensor_elem(nd->left) + _sum_of_tensor_elem(nd->right);
    }

  template <typename uni10_type>
    size_t Network<uni10_type>::memory_requirement(){
      if(rollcall() >= 0)
        return 0;
      size_t usage = 0;
      for(int i = 0; i < (int)leafs.size(); i++)
        usage += leafs[i]->elemNum;
      usage *= 2;
      size_t max_usage = 0;
      _elem_usage(root, usage, max_usage);
      return max_usage * sizeof(uni10_type);
    }

  template <typename uni10_type>
    size_t Network<uni10_type>::_elem_usage(Node<uni10_type>* nd, size_t& usage, size_t& max_usage)const{
      if(nd == NULL)
        return 0;
      size_t child_usage = _elem_usage(nd->left, usage, max_usage) + _elem_usage(nd->right, usage, max_usage);
      usage += nd->elemNum;
      max_usage = usage > max_usage ? usage : max_usage;
      usage -= child_usage;
      return nd->elemNum;
    }

  template <typename uni10_type>
    size_t Network<uni10_type>::max_tensor_elemNum(){
      if(rollcall() >= 0)
        return 0;
      size_t max_num = 0;
      Node<uni10_type> max_nd;
      _max_tensor_elemNum(root, max_num, max_nd);
      return max_num;
    }

  template <typename uni10_type>
    void Network<uni10_type>::_max_tensor_elemNum(Node<uni10_type>* nd, size_t& max_num, Node<uni10_type>& max_nd) const{
      if(nd == NULL)
        return;
      _max_tensor_elemNum(nd->left, max_num, max_nd);
      _max_tensor_elemNum(nd->right, max_num, max_nd);
      if(nd->elemNum > max_num){
        max_num = nd->elemNum;
        max_nd = *nd;
      }
    }

  template <typename uni10_type>
    std::string Network<uni10_type>::profile(bool print){
      std::ostringstream os;
      int miss;
      if((miss = rollcall()) >= 0){
        os<<"\nTensor '"<<names[miss]<<"' has not yet been given!\n\n";
        if(print){
          std::cout<<os.str();
          return "";
        }
        return os.str();
      }
      os<<"\n===== Network profile =====\n";
      os<<"Memory Requirement: "<<memory_requirement()<<std::endl;
      //os<<"Sum of memory usage: "<<sum_of_memory_usage()<<std::endl;
      size_t max_num = 0;
      Node<uni10_type> max_nd;
      _max_tensor_elemNum(root, max_num, max_nd);
      os<<"Maximun tensor: \n";
      os<<"  elemNum: "<<max_num<<"\n  "<<max_nd.labels.size()<<" bonds and labels: ";
      for(int i = 0; i < (int)max_nd.labels.size(); i++)
        os<< max_nd.labels[i] << ", ";
      os<<std::endl;
      os<<"===========================\n\n";
      if(print){
        std::cout<<os.str();
        return "";
      }
      return os.str();
    }

  template <typename uni10_type>
    void Network<uni10_type>::preprint(std::ostream& os, Node<uni10_type>* nd, int layer)const{
      if(nd == NULL)
        return;
      for(int i = 0; i < layer; i++)
        os<<"|   ";
      if(nd->T)
        os<<nd->name << "(" << nd->elemNum << "): ";
      else
        os<<"*("<<nd->elemNum<<"): ";
      for(int i = 0; i < (int)nd->labels.size(); i++)
        os<< nd->labels[i] << ", ";
      os<<std::endl;
      preprint(os, nd->left, layer+1);
      preprint(os, nd->right, layer+1);
    }

  template <typename uni10_type>
    void Network<uni10_type>::findConOrd(Node<uni10_type>* nd){
      if(nd == NULL || conOrder.size() == tensors.size())
        return;
      if(nd->T){
        bool found = false;
        for(int i = 0; i < (int)tensors.size(); i++)
          if(nd->T == tensors[i]){
            conOrder.push_back(i);
            found = true;
            break;
          }
        uni10_error_msg(!found, "%s","Fatal error(code = N2). Please contact the developer of the uni10 library.");
      }
      findConOrd(nd->left);
      findConOrd(nd->right);
    }

  template <typename uni10_type>
    void Network<uni10_type>::addSwap(){
      int Tnum = leafs.size();
      findConOrd(root);
      uni10_error_msg(!(Tnum == (int)conOrder.size()), "%s", "Fatal error(code = N3). Please contact the developer of the uni10 library.");
      std::vector<int> tenOrder = conOrder;

      /*
      std::cout << "conOrder: ";

      for(uni10_uint64 i = 0; i < tenOrder.size(); i++)
        std::cout << tenOrder[i] << " ";
      std::cout << std::endl;

      for(uni10_uint64 i = 0; i < tenOrder.size(); i++)
        std::cout << names[tenOrder[i]] << " ";
      std::cout << std::endl;
      */

      std::vector<_Swap> tenSwaps = recSwap(tenOrder);
      std::vector<_Swap> swtmp;
      for(int s = 0; s < (int)tenSwaps.size(); s++){
        swtmp = tensors[tenSwaps[s].b1]->exSwap(*(tensors[tenSwaps[s].b2]));
        swaps_arr[tenSwaps[s].b1].insert(swaps_arr[tenSwaps[s].b1].end(), swtmp.begin(), swtmp.end());
      }
      //Distinct Swaps of each tensors
      for(int t = 0; t < Tnum; t++){
        std::map<int, bool> recs;
        std::map<int, bool>::iterator it;
        int bondNum = (int)tensors[t]->bonds->size();
        int is, ib;
        //int hash;
        for(int s = 0; s < (int)swaps_arr[t].size(); s++){
          if(swaps_arr[t][s].b1 < swaps_arr[t][s].b2){
            is = swaps_arr[t][s].b1;
            ib = swaps_arr[t][s].b2;
          }
          else{
            is = swaps_arr[t][s].b2;
            ib = swaps_arr[t][s].b1;
          }
          int hash = is * bondNum + ib;
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
    }

  template class Network<uni10_double64>;
  template class Network<uni10_complex128>;

}; /* namespace uni10 */
