#include <limits.h>

#include "uni10/uni10_api/network_tools/PseudoTensor.h"

namespace uni10{

  PseudoTensor::PseudoTensor(): bit(0), cost(0), is_new(true), max_label(INT_MIN){};

  PseudoTensor::PseudoTensor(std::vector<std::string>& _order, std::map<uni10_int32, uni10_int32>& _label_dim, uni10_uint64 _bit, uni10_float32 _cost, uni10_bool _is_new):
  order(_order), bit(_bit), cost(_cost), is_new(_is_new), label_dim(_label_dim){

    max_label = label_dim.rbegin()->first;

  }

  void PseudoTensor::printPseudoTensor() const {

    fprintf(stdout, "\n");
    fprintf(stdout, "=======================\n");
    for(uni10_int32 i = 0; i < (uni10_int32)order.size(); i++){
      if(i == 0)
        fprintf(stdout, "ORDER: %s", order[i].c_str());
      else
        fprintf(stdout, ", %s", order[i].c_str());
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "Bit : %ld\n", bit);
    fprintf(stdout, "Cost: %.7f\n", cost);
    if(is_new)
      fprintf(stdout, "New : True\n");
    else
      fprintf(stdout, "New : False\n");

    if(label_dim.size() != 0){

      fprintf(stdout, "{label, dim}: ");
      std::map<uni10_int32, uni10_int32>::const_iterator it = label_dim.begin();
      fprintf(stdout, "{%d, %d}", it->first, it->second);
      for(it++; it != label_dim.end(); it++)
        fprintf(stdout, ", {%d, %d}", it->first, it->second);
      fprintf(stdout, "\n");
      fprintf(stdout, "Max_label: %d\n", max_label);

    }else
      fprintf(stdout, "Scaler !!\n");
    fprintf(stdout, "=======================\n");
    fprintf(stdout, "\n");

  }

  void PseudoTensor::netorder(char** netorder){

    uni10_int32 num_of_sub_tree;
    uni10_int32* size_of_sub_tree = NULL;
    char* buf_netorder = NULL;
    char** Tnames = NULL;
    uni10_int32 numT;

    uni10_int32 netorder_len = 0;
    uni10_int32 unit_ptr_size = 128;
    uni10_int32 netorder_size = unit_ptr_size;
    uni10_int32 cnt_ten = 0;
    buf_netorder = (char*)calloc(unit_ptr_size, sizeof(char));

    bufferInfo(&num_of_sub_tree, &size_of_sub_tree, &Tnames, &numT);

    for(uni10_int32 g = 0; g < (uni10_int32)num_of_sub_tree; g++){
      if( g == 0 ){
        netorder_len += num_of_sub_tree-1;
        if(netorder_len > netorder_size){
          netorder_size = (netorder_len/unit_ptr_size+1)*unit_ptr_size;
          buf_netorder = (char*)realloc(buf_netorder, netorder_size*sizeof(char));
        }
        for(uni10_int32 n = 0; n < num_of_sub_tree-1; n++)
          strcat(buf_netorder, "(");
      }

      for(uni10_int32 t = 0; t < size_of_sub_tree[g]; t++){
        if(t == 0){
          netorder_len += size_of_sub_tree[g]-1 + strlen(Tnames[cnt_ten]);
          if(netorder_len > netorder_size){
            netorder_size = (netorder_len/unit_ptr_size+1)*unit_ptr_size;
            buf_netorder = (char*)realloc(buf_netorder, netorder_size*sizeof(char));
          }

          for(uni10_int32 n = 0; n < size_of_sub_tree[g]-1; n++)
            strcat(buf_netorder, "(");

          strcat(buf_netorder, Tnames[cnt_ten]);
          cnt_ten++;

        }else{
          netorder_len += 2 + strlen(Tnames[cnt_ten]);
          if(netorder_len > netorder_size){
            netorder_size = (netorder_len/unit_ptr_size+1)*unit_ptr_size;
            buf_netorder = (char*)realloc(buf_netorder, netorder_size*sizeof(char));
          }

          strcat(buf_netorder, " ");
          strcat(buf_netorder, Tnames[cnt_ten]);
          strcat(buf_netorder, ")");
          cnt_ten++;

        }
      }

      if(g != 0){
        netorder_len += 1;
        if(netorder_len > netorder_size){
          netorder_size = (netorder_len/unit_ptr_size+1)*unit_ptr_size;
          buf_netorder = (char*)realloc(buf_netorder, netorder_size*sizeof(char));
        }
        strcat(buf_netorder, ")");
      }

    }

    *netorder = buf_netorder;

    if(size_of_sub_tree != NULL)
      free(size_of_sub_tree);
    if(Tnames != NULL){
      for(uni10_int32 i = 0; i < numT; i++)
        free(Tnames[i]);
      free(Tnames);
    }

  }

  void PseudoTensor::bufferInfo(uni10_int32* num_of_sub_tree, uni10_int32** size_of_sub_tree, char*** Tnames, uni10_int32 *numT){

    uni10_int32 unit_int_ptr_size = 12;
    uni10_int32 unit_double_cptr_size = 32;
    uni10_int32 int_ptr_size = unit_int_ptr_size;
    uni10_int32 double_cptr_size = unit_double_cptr_size;

    uni10_int32* buf_sub_tree_size = (uni10_int32*)calloc(int_ptr_size, sizeof(uni10_int32));
    char** buf_Tnames = (char**)calloc(double_cptr_size, sizeof(char*));
    *num_of_sub_tree = 0;
    uni10_int32 cnt_ten = 0;

    for(uni10_int32 i = order.size()-1; i >= 0; --i){

      if( i == (uni10_int32)order.size()-1 || (order[i] == "#" && order[i+1] != "#")){
        *num_of_sub_tree += 1;
        if(*num_of_sub_tree > int_ptr_size){
          int_ptr_size = (*num_of_sub_tree/unit_int_ptr_size+1)*unit_int_ptr_size;
          buf_sub_tree_size = (uni10_int32*)realloc(buf_sub_tree_size, int_ptr_size*sizeof(uni10_int32));
        }
      }

      else if( order[i] != "#"){
        buf_sub_tree_size[*num_of_sub_tree-1]++;
        cnt_ten++;
        if(cnt_ten > double_cptr_size){
          double_cptr_size = (cnt_ten/unit_double_cptr_size+1)*unit_double_cptr_size;
          buf_Tnames = (char**)realloc(buf_Tnames, double_cptr_size*sizeof(char*));
        }
        buf_Tnames[cnt_ten-1] = (char*)malloc(order[i].size()*sizeof(char));
        strcpy(buf_Tnames[cnt_ten-1], order[i].c_str());
      }
    }

    *size_of_sub_tree = buf_sub_tree_size;
    *Tnames = buf_Tnames;
    *numT = cnt_ten;

  }


};
