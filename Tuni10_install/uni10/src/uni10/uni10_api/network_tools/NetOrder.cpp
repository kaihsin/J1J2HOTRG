#include <float.h>

#include "uni10/uni10_api/network_tools/NetOrder.h"

namespace uni10{

  NetOrder::NetOrder(): tensor_set(ary2d_pten()), xi_min(-1.), numT(0), netorder(NULL){}

  NetOrder::NetOrder(const ary1d_uptr_d& tens): xi_min(-1.), netorder(NULL){

    tensor_set.assign(tens.size(), ary1d_pten());
    numT = tens.size();

    uni10_int32 cnt = 0;
    for(uni10_int32 t = 0; t < (uni10_int32)tens.size(); t++){

      std::vector<uni10_int32> label = tens[t]->label();

      PseudoTensor tmp;
      tmp.order.push_back(tens[t]->getName());

      tmp.bit = 0;
      tmp.bit += pow(2, cnt);
      cnt++;

      tmp.cost = 0.0;

      for(uni10_int32 i = 0; i < (uni10_int32)label.size(); i++){
        uni10_int32 BDim = tens[t]->bond(i).dim();
        if(xi_min < 0 || xi_min > BDim)
          xi_min = BDim;
        tmp.label_dim[label[i]] = BDim;

      }

      tmp.max_label = tmp.label_dim.rbegin()->first;
      tensor_set[0].push_back(tmp);
    }

  }

  NetOrder::NetOrder(const ary1d_uptr_z& tens): xi_min(-1.), netorder(NULL){

    tensor_set.assign(tens.size(), ary1d_pten());
    numT = tens.size();

    uni10_int32 cnt = 0;
    for(uni10_int32 t = 0; t < (uni10_int32)tens.size(); t++){

      std::vector<uni10_int32> label = tens[t]->label();

      PseudoTensor tmp;
      tmp.order.push_back(tens[t]->getName());

      tmp.bit = 0;
      tmp.bit += pow(2, cnt);
      cnt++;

      tmp.cost = 0.0;

      for(uni10_int32 i = 0; i < (uni10_int32)label.size(); i++){
        uni10_int32 BDim = tens[t]->bond(i).dim();
        if(xi_min < 0 || xi_min > BDim)
          xi_min = BDim;
        tmp.label_dim[label[i]] = BDim;

      }

      tmp.max_label = tmp.label_dim.rbegin()->first;
      tensor_set[0].push_back(tmp);
    }

  }

  NetOrder::~NetOrder(){

    if(netorder != NULL)
      free(netorder);

  }

  char* NetOrder::generate_order(){

    uni10_float32 mu_cap = 1.0;
    uni10_float32 mu_old = 0.0;
    uni10_float32 mu_next;

    while(tensor_set.back().size() == 0){

      mu_next = FLT_MAX;

      for(int c = 1; c < numT; c++){

        for(int d1 = 0; d1 < (c+1)/2; d1++){
          uni10_int32 d2 = c -d1 -1;
          uni10_int32 n1 = tensor_set[d1].size();
          uni10_int32 n2 = tensor_set[d2].size();
          for(int i1 = 0; i1 < n1; i1++){
            uni10_int32 i2_start = d1==d2 ? i1+1 : 0;
            for(int i2 = i2_start; i2 < n2; i2++){
              PseudoTensor t1 = tensor_set[d1][i1];
              PseudoTensor t2 = tensor_set[d2][i2];
              if(this->is_disjoint(t1, t2))
                continue;
              if(this->is_overlap(t1, t2))
                continue;
              uni10_float32 mu = get_cost(t1, t2);
              uni10_float32 mu_0 = (t1.is_new || t2.is_new) ? 0.0 : mu_old;

              if(mu > mu_cap && mu < mu_next)
                mu_next = mu;
              if(mu > mu_0 && mu <= mu_cap){
                PseudoTensor t3 = psesudocontract(t1, t2);
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
      for(uni10_int32 s = 0; s < (uni10_int32)tensor_set.size(); s++)
        for(uni10_int32 t = 0; t < (uni10_int32)tensor_set[s].size(); t++)
          tensor_set[s][t].is_new = false;

    }

    tensor_set.back()[0].netorder(&netorder);
    //for(int i = 0; i < (int)tensor_set.back()[0].order.size(); i++)
    //  printf("%s ", tensor_set.back()[0].order[i].c_str());
    //printf("\n");

    //exit(0);

    return netorder;

  }

  uni10_bool NetOrder::is_disjoint(const PseudoTensor& t1, const PseudoTensor& t2){

    uni10_bool isdisjoint = true;

    std::map<uni10_int32, uni10_int32>::const_iterator it1 = t1.label_dim.begin();
    std::map<uni10_int32, uni10_int32>::const_iterator it2 = t2.label_dim.begin();

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

  uni10_bool NetOrder::is_overlap(const PseudoTensor& t1, const PseudoTensor& t2){

    return (t1.bit & t2.bit) > 0;

  }

  PseudoTensor NetOrder::psesudocontract(const PseudoTensor& t1, const PseudoTensor& t2){

    PseudoTensor t3;
    assert(!is_disjoint(t1, t2));

    t3.order.reserve(t1.order.size()+t2.order.size());
    t3.order.insert(t3.order.end(), t1.order.begin(), t1.order.end());
    t3.order.insert(t3.order.end(), t2.order.begin(), t2.order.end());
    t3.order.push_back("#");

    std::map<uni10_int32, uni10_int32>::const_iterator it1  = t1.label_dim.begin();
    std::map<uni10_int32, uni10_int32>::const_iterator it2  = t2.label_dim.begin();

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

    std::map<uni10_int32, uni10_int32>::const_iterator itMax = t1.max_label > t2.max_label ? it1 : it2;
    std::map<uni10_int32, uni10_int32>::const_iterator itMaxE = t1.max_label > t2.max_label ? t1.label_dim.end() : t2.label_dim.end();

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

  uni10_float32 NetOrder::get_cost(const PseudoTensor& t1, const PseudoTensor& t2){

    uni10_float32 cost = 1.;

    std::map<uni10_int32, uni10_int32>::const_iterator it1 = t1.label_dim.begin();
    std::map<uni10_int32, uni10_int32>::const_iterator it2 = t2.label_dim.begin();

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

    std::map<int, int>::const_iterator itMax  = t1.max_label > t2.max_label ? it1 : it2;
    std::map<int, int>::const_iterator itMaxE = t1.max_label > t2.max_label ? t1.label_dim.end() : t2.label_dim.end();

    while(itMax != itMaxE){
      cost *= itMax->second;	
      itMax++;
    }

    cost += t1.cost + t2.cost;

    return cost;
  }

}
