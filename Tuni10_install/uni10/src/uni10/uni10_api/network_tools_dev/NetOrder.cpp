#include <float.h>

#include "uni10/uni10_api/network_tools_dev/NetOrder.h"

namespace uni10{

  NetOrder_dev::NetOrder_dev(): tensor_set(ary2d_pten()), xi_min(-1.), numT(0), names(NULL){}

  NetOrder_dev::~NetOrder_dev(){

  }

  NetOrder_dev::NetOrder_dev(const ary1d_uptr& tens_type, const ary2d_label& label_arr, const ary1d_name& _names): xi_min(-1.), names(& _names){

    tensor_set.assign(tens_type.size(), ary1d_pten());
    numT = tens_type.size();

    uni10_int cnt = 0;
    uni10_int typeID;

    for(uni10_int t = 0; t < (uni10_int)tens_type.size(); t++){

      PseudoTensor_dev tmp;
      typeID = tens_type[t].second;

      // The order of tensor idx in the tensor list.
      tmp.order_idx.push_back(t);

      tmp.bit = 0;
      tmp.bit += pow(2, cnt);
      cnt++;

      tmp.cost = 0.0;

      // Counting the lable_dim.
      for(uni10_int i = 0; i < (uni10_int)label_arr[t].size(); i++){

        uni10_int BDim = (typeID == 1) ? 
          ((UniTensor<uni10_double64>*)tens_type[t].first)->bond(i).dim() : ((UniTensor<uni10_complex128>*)tens_type[t].first)->bond(i).dim();

        if(xi_min < 0 || xi_min > BDim)
          xi_min = BDim;
        tmp.label_dim[ label_arr[t][i] ] = BDim;

      }

      tmp.max_label = tmp.label_dim.rbegin()->first;

      tensor_set[0].push_back(tmp);

    }

  }

  ary1d_order NetOrder_dev::generate_order(){

    uni10_float32 mu_cap = 1.0;
    uni10_float32 mu_old = 0.0;
    uni10_float32 mu_next;

    while(tensor_set.back().size() == 0){

      mu_next = FLT_MAX;

      for(int c = 1; c < numT; c++){

        for(int d1 = 0; d1 < (c+1)/2; d1++){
          uni10_int d2 = c -d1 -1;
          uni10_int n1 = tensor_set[d1].size();
          uni10_int n2 = tensor_set[d2].size();
          for(int i1 = 0; i1 < n1; i1++){
            uni10_int i2_start = d1==d2 ? i1+1 : 0;
            for(int i2 = i2_start; i2 < n2; i2++){
              PseudoTensor_dev t1 = tensor_set[d1][i1];
              PseudoTensor_dev t2 = tensor_set[d2][i2];
              if(this->is_disjoint(t1, t2))
                continue;
              if(this->is_overlap(t1, t2))
                continue;
              uni10_float32 mu = get_cost(t1, t2);
              uni10_float32 mu_0 = (t1.is_new || t2.is_new) ? 0.0 : mu_old;

              if(mu > mu_cap && mu < mu_next)
                mu_next = mu;
              if(mu > mu_0 && mu <= mu_cap){
                PseudoTensor_dev t3 = psesudocontract(t1, t2);
                uni10_bool exsist = false;
                for(int i = 0; i < (int)tensor_set[c].size(); i++){
                  if(t3.bit == tensor_set[c][i].bit){
                    if(t3.cost < tensor_set[c][i].cost)
                      tensor_set[c][i] = t3;
                    exsist = true;
                    break;
                  }

                }

                if(!exsist)
                  tensor_set[c].push_back(t3);

              }

            }

          }

        }

      }

      mu_old = mu_cap;
      mu_cap = std::max(mu_next, mu_cap*xi_min);
      for(uni10_int s = 0; s < (uni10_int)tensor_set.size(); s++)
        for(uni10_int t = 0; t < (uni10_int)tensor_set[s].size(); t++)
          tensor_set[s][t].is_new = false;

    }

    return tensor_set.back()[0].order_idx;

  }

  uni10_bool NetOrder_dev::is_disjoint(const PseudoTensor_dev& t1, const PseudoTensor_dev& t2){

    uni10_bool isdisjoint = true;

    std::map<uni10_int, uni10_int>::const_iterator it1 = t1.label_dim.begin();
    std::map<uni10_int, uni10_int>::const_iterator it2 = t2.label_dim.begin();

    while(it1 != t1.label_dim.end() && it2 != t2.label_dim.end()){

      if(it1->first < it2->first)
        ++it1;
      else if(it2->first < it1->first)
        ++it2;
      else{
        isdisjoint = false;
        break;	
      }

    }

    return isdisjoint;

  }

  uni10_bool NetOrder_dev::is_overlap(const PseudoTensor_dev& t1, const PseudoTensor_dev& t2){

    return (t1.bit & t2.bit) > 0;

  }

  PseudoTensor_dev NetOrder_dev::psesudocontract(const PseudoTensor_dev& t1, const PseudoTensor_dev& t2){

    PseudoTensor_dev t3;
    assert(!is_disjoint(t1, t2));

    //t3.order_str.reserve(t1.order_str.size()+t2.order_str.size());
    //t3.order_str.insert(t3.order_str.end(), t1.order_str.begin(), t1.order_str.end());
    //t3.order_str.insert(t3.order_str.end(), t2.order_str.begin(), t2.order_str.end());
    //t3.order_str.push_back("#");

    t3.order_idx.reserve(t1.order_idx.size()+t2.order_idx.size());
    t3.order_idx.insert(t3.order_idx.end(), t1.order_idx.begin(), t1.order_idx.end());
    t3.order_idx.insert(t3.order_idx.end(), t2.order_idx.begin(), t2.order_idx.end());
    t3.order_idx.push_back(-1);

    std::map<uni10_int, uni10_int>::const_iterator it1  = t1.label_dim.begin();
    std::map<uni10_int, uni10_int>::const_iterator it2  = t2.label_dim.begin();

    while(it1 != t1.label_dim.end() && it2 != t2.label_dim.end()){

      if(it1->first < it2->first){
        t3.label_dim[it1->first] = it1->second;
        ++it1;
      }
      else if(it1->first > it2->first){
        t3.label_dim[it2->first] = it2->second;
        ++it2;
      }
      else{
        ++it1;
        ++it2;
      }

    }

    std::map<uni10_int, uni10_int>::const_iterator itMax  = t1.max_label > t2.max_label ? it1 : it2;
    std::map<uni10_int, uni10_int>::const_iterator itMaxE = t1.max_label > t2.max_label ? t1.label_dim.end() : t2.label_dim.end();

    while(itMax != itMaxE){
      t3.label_dim[itMax->first] = itMax->second;
      itMax++;
    }

    t3.cost = get_cost(t1, t2);
    t3.bit = t1.bit ^ t2.bit;
    t3.is_new = true;

    if(t3.label_dim.size() != 0)
      t3.max_label = t3.label_dim.rbegin()->first;

    return t3;

  }

  uni10_float32 NetOrder_dev::get_cost(const PseudoTensor_dev& t1, const PseudoTensor_dev& t2){

    uni10_float32 cost = 1.;

    std::map<uni10_int, uni10_int>::const_iterator it1 = t1.label_dim.begin();
    std::map<uni10_int, uni10_int>::const_iterator it2 = t2.label_dim.begin();

    while(it1 != t1.label_dim.end() && it2 != t2.label_dim.end()){

      if(it1->first < it2->first){
        cost *= it1->second;
        ++it1;
      }
      else if(it2->first < it1->first){
        cost *= it2->second;
        ++it2;
      }
      else{
        cost *= it1->second;
        ++it1;
        ++it2;
      }

    }

    std::map<uni10_int, uni10_int>::const_iterator itMax  = t1.max_label > t2.max_label ? it1 : it2;
    std::map<uni10_int, uni10_int>::const_iterator itMaxE = t1.max_label > t2.max_label ? t1.label_dim.end() : t2.label_dim.end();

    while(itMax != itMaxE){
      cost *= itMax->second;	
      itMax++;
    }

    cost += t1.cost + t2.cost;

    return cost;
  }

}
