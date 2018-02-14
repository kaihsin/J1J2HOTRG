/****************************************************************************
*  @file Bond.cpp
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
*  @brief Implementation file for Matrix class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_api/Bond.h"

namespace uni10{

  Bond::Bond(bondType _type, size_t dim) : m_type(_type){
    Qnum q0(0);
    std::vector<Qnum> qnums(dim, q0);
    setting(qnums);
  }

  Bond::Bond(bondType _type, const std::vector<Qnum>& qnums) : m_type(_type){
    setting(qnums);
  }

  Bond::Bond(const Bond& _b):m_type(_b.m_type), m_dim(_b.m_dim), Qnums(_b.Qnums), Qdegs(_b.Qdegs), offsets(_b.offsets){
  }

  bondType Bond::type()const{
    return m_type;
  }
  uni10_int32 Bond::dim()const{
    return m_dim;
  }

  void Bond::assign(bondType _type, size_t dim){

    m_type = _type;
    Qnums.clear();
    Qdegs.clear();
    offsets.clear();
    Qnum q0(0);
    std::vector<Qnum> qnums(dim, q0);
    setting(qnums);

  }

  void Bond::assign(bondType _type, const std::vector<Qnum>& qnums){

    m_type = _type;
    Qnums.clear();
    Qdegs.clear();
    offsets.clear();
    setting(qnums);

  }

  void Bond::setting(const std::vector<Qnum>& qnums){

    uni10_error_msg(!(qnums.size() > 0), "%s", "Cannot create a bond of dimension 0.");

    std::map<Qnum, bool> mark;

    uni10_int32 cnt = 0;

    m_dim = 0;

    for(uni10_int32 i = 0; i < (uni10_int32)qnums.size(); i++){
      if(i == 0 || !(qnums[i] == qnums[i - 1])){
        Qnums.push_back(qnums[i]);
        Qdegs.push_back(1);
        offsets.push_back(m_dim);
        cnt++;
      }
      else
        Qdegs[cnt - 1]++;
      m_dim++;
    }
  }

  Bond::~Bond(){}

  std::ostream& operator<< (std::ostream& os, const Bond& b){
    if(b.m_type == BD_IN)
      os<<"IN : ";
    else
      os<<"OUT: ";
    for(uni10_int32 i = 0; i < (uni10_int32)b.Qnums.size(); i++)
      os << b.Qnums[i] << "|" << b.Qdegs[i]<<", ";
    os<<"Dim = "<< b.m_dim << std::endl;
    return os;
  }

  std::map<Qnum, uni10_int32> Bond::degeneracy()const{
    std::map<Qnum, uni10_int32>hst;
    for(uni10_int32 i = 0; i < (uni10_int32)Qnums.size(); i++){
      if(hst.find(Qnums[i]) == hst.end())
        hst[Qnums[i]] = Qdegs[i];
      else
        hst[Qnums[i]] += Qdegs[i];
    }
    return hst;
  }

  std::vector<Qnum> Bond::Qlist()const{
    std::vector<Qnum>list(m_dim);
    uni10_int32 cnt = 0;
    for(uni10_int32 q = 0; q < (uni10_int32)Qnums.size(); q++)
      for(uni10_int32 d = 0; d < Qdegs[q]; d++){
        list[cnt] = Qnums[q];
        cnt++;
      }
    return list;
  }

  bool operator== (const Bond& b1, const Bond& b2){
    return (b1.m_type == b2.m_type) && (b1.Qnums == b2.Qnums) && (b1.Qdegs == b2.Qdegs);
  }
  Bond& Bond::change(bondType tp){
    if(m_type != tp){
      for(uni10_int32 q = 0; q < (uni10_int32)Qnums.size(); q++)
        Qnums[q] = -Qnums[q];
      m_type = tp;
    }
    return *this;
  }
  Bond& Bond::dummy_change(bondType tp){
    if(m_type != tp)
      m_type = tp;
    return *this;
  }
  Bond& Bond::combine(Bond bd){

    bd.change(m_type);
    std::vector<Qnum> qnums;
    std::vector<uni10_int32> qdegs;
    offsets.clear();
    m_dim = 0;
    Qnum qnum;
    uni10_int32 qdim;
    uni10_int32 cnt = 0;
    for(uni10_int32 q = 0; q < (uni10_int32)Qnums.size(); q++)
      for(uni10_int32 d = 0; d < Qdegs[q]; d++){
        for(uni10_int32 qq = 0; qq < (uni10_int32)bd.Qnums.size(); qq++){
          qnum = Qnums[q] * bd.Qnums[qq];
          qdim = bd.Qdegs[qq];
          if(qnums.size() == 0 || !(qnum == qnums[cnt - 1])){
            qnums.push_back(qnum);
            qdegs.push_back(qdim);
            offsets.push_back(m_dim);
            cnt++;
          }
          else
            qdegs[cnt - 1] += qdim;
          m_dim += qdim;
        }
      }
    Qnums = qnums;
    Qdegs = qdegs;

    return *this;
  }

  Bond combine(bondType tp, const std::vector<Bond>& bds){

    uni10_error_msg((bds.size() == 0), "%s", "There should be at least one bond in the input vector to be combined.");

    if(bds.size() == 1){
      Bond bd = bds[0];
      bd.change(tp);
      return bd;
    }
    uni10_int32 bd_num = bds.size();
    Bond outBond1 = bds[bd_num - 1];
    Bond outBond2 = bds[bd_num - 2];
    uni10_int32 b = 0;
    outBond2.change(tp);
    outBond2.combine(outBond1);
    for(b = 0; b < bd_num - 2; b++){
      if(b % 2 == 0){
        outBond1 = bds[bd_num - 3 - b];
        outBond1.change(tp);
        outBond1.combine(outBond2);
      }
      else{
        outBond2 = bds[bd_num - 3 - b];
        outBond2.change(tp);
        outBond2.combine(outBond1);
      }
    }
    if(b % 2 == 0)
      return outBond2;
    else
      return outBond1;

  }

  Bond combine(const std::vector<Bond>& bds){
    uni10_error_msg(bds.size() == 0, "%s", "There should be at least one bond in the input vector to be combined.");
    return combine(bds[0].m_type, bds);
  }

};	/* namespace uni10 */
