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
#include <algorithm>

#include "uni10/uni10_api/network_tools/Node.h"

namespace uni10{

  template <typename uni10_type>
    Node<uni10_type>::Node(): T(NULL), elemNum(0), parent(NULL), left(NULL), right(NULL), point(0){

    }

  template <typename uni10_type>
    Node<uni10_type>::Node(UniTensor<uni10_type>* Tp): T(Tp), labels(*Tp->labels), bonds(*Tp->bonds), elemNum(*Tp->U_elemNum), name(*Tp->name), parent(NULL), left(NULL), right(NULL), point(0){

      uni10_error_msg(!(*Tp->status & Tp->HAVEBOND), "%s", "Cannot create node of a network from tensor without bond.");

    }

  template <typename uni10_type>
    Node<uni10_type>::Node(const Node<uni10_type>& nd): T(nd.T), labels(nd.labels), bonds(nd.bonds), elemNum(nd.elemNum), parent(nd.parent), left(nd.left), right(nd.right), point(nd.point){
    }

  template <typename uni10_type>
    Node<uni10_type>::Node(std::vector<Bond>& _bonds, std::vector<uni10_int32>& _labels): T(NULL), labels(_labels), bonds(_bonds), parent(NULL), left(NULL), right(NULL), point(0){
      elemNum = cal_elemNum(bonds);
    }

  template <typename uni10_type>
    Node<uni10_type>::~Node(){

    }

  template <typename uni10_type>
    void Node<uni10_type>::delink(){
      parent = NULL;
      left = NULL;
      right = NULL;
      point = 0;
    }

  template <typename uni10_type>
    Node<uni10_type> Node<uni10_type>::contract(Node<uni10_type>* nd){
      int AbondNum = bonds.size();
      int BbondNum = nd->bonds.size();
      std::vector<Bond> cBonds;
      std::vector<int> markB(BbondNum, 0);
      std::vector<int> newLabelC;
      int conBondNum = 0;
      bool match;
      for(int a = 0; a < AbondNum; a++){
        match = false;
        for(int b = 0; b < BbondNum; b++)
          if(labels[a] == nd->labels[b]){
            markB[b] = 1;
            match = true;
            conBondNum++;
            break;
          }
        if(!match){
          newLabelC.push_back(labels[a]);
          cBonds.push_back(bonds[a]);
        }
      }
      for(int b = 0; b < BbondNum; b++)
        if(markB[b] == 0){
          newLabelC.push_back(nd->labels[b]);
          cBonds.push_back(nd->bonds[b]);
        }
      int rBondNum = AbondNum - conBondNum;
      int cBondNum = BbondNum - conBondNum;
      for(int a = 0; a < rBondNum; a++)
        cBonds[a].change(BD_IN);
      for(int a = 0; a < cBondNum; a++)
        cBonds[rBondNum + a].change(BD_OUT);

      //Node par(cBonds, newLabelC);
      return Node(cBonds, newLabelC);
    }

  template<typename uni10_type>
    float Node<uni10_type>::metric(Node<uni10_type>* nd){	//Bigger is better
      int AbondNum = bonds.size();
      int BbondNum = nd->bonds.size();
      std::vector<Bond> cBonds;
      std::vector<int> markB(BbondNum, 0);
      int conBondNum = 0;
      bool match;
      for(int a = 0; a < AbondNum; a++){
        match = false;
        for(int b = 0; b < BbondNum; b++)
          if(labels[a] == nd->labels[b]){
            markB[b] = 1;
            match = true;
            conBondNum++;
            break;
          }
        if(!match)
          cBonds.push_back(bonds[a]);
      }
      if(conBondNum == 0)
        return -1;
      for(int b = 0; b < BbondNum; b++)
        if(markB[b] == 0)
          cBonds.push_back(nd->bonds[b]);
      int rBondNum = AbondNum - conBondNum;
      int cBondNum = BbondNum - conBondNum;
      for(int a = 0; a < rBondNum; a++)
        cBonds[a].change(BD_IN);
      for(int a = 0; a < cBondNum; a++)
        cBonds[rBondNum + a].change(BD_OUT);
      int64_t newElemNum = cal_elemNum(cBonds);
      return float(elemNum + nd->elemNum) / newElemNum;
    }

  template <typename uni10_type>
    uni10_uint64 Node<uni10_type>::cal_elemNum(std::vector<Bond>& _bonds){
      int rBondNum = 0;
      int cBondNum = 0;
      for(int b = 0; b < (int)_bonds.size(); b++)
        if(_bonds[b].type() == BD_IN)
          rBondNum++;
        else if(_bonds[b].type() == BD_OUT)
          cBondNum++;
      Qnum qnum(0, PRT_EVEN);
      size_t dim;
      std::map<Qnum,size_t> row_QnumMdim;
      std::vector<int> row_offs(rBondNum, 0);
      if(rBondNum){
        while(1){
          qnum.assign();
          dim = 1;
          for(int b = 0; b < rBondNum; b++){
            qnum = qnum * _bonds[b].const_getQnums()[row_offs[b]];
            dim *= _bonds[b].const_getQdegs()[row_offs[b]];
          }
          if(row_QnumMdim.find(qnum) != row_QnumMdim.end())
            row_QnumMdim[qnum] += dim;
          else
            row_QnumMdim[qnum] = dim;
          int bidx;
          for(bidx = rBondNum - 1; bidx >= 0; bidx--){
            row_offs[bidx]++;
            if(row_offs[bidx] < (int)_bonds[bidx].const_getQnums().size())
              break;
            else
              row_offs[bidx] = 0;
          }
          if(bidx < 0)	//run over all row_bond offsets
            break;
        }
      }
      else{
        qnum.assign();
        row_QnumMdim[qnum] = 1;
      }

      std::map<Qnum,size_t> col_QnumMdim;
      std::vector<int> col_offs(cBondNum, 0);
      if(cBondNum){
        while(1){
          qnum.assign();
          dim = 1;
          for(int b = 0; b < cBondNum; b++){
            qnum = qnum * _bonds[b + rBondNum].const_getQnums()[col_offs[b]];
            dim *= _bonds[b + rBondNum].const_getQdegs()[col_offs[b]];
          }
          if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
            if(col_QnumMdim.find(qnum) != col_QnumMdim.end())
              col_QnumMdim[qnum] += dim;
            else
              col_QnumMdim[qnum] = dim;
          }
          int bidx;
          for(bidx = cBondNum - 1; bidx >= 0; bidx--){
            col_offs[bidx]++;
            if(col_offs[bidx] < _bonds[bidx + rBondNum].const_getQnums().size())
              break;
            else
              col_offs[bidx] = 0;
          }
          if(bidx < 0)	//run over all row_bond offsets
            break;
        }
      }
      else{
        qnum.assign();
        if(row_QnumMdim.find(qnum) != row_QnumMdim.end())
          col_QnumMdim[qnum] = 1;
      }
      size_t _elemNum = 0;
      std::map<Qnum,size_t>::iterator it;
      std::map<Qnum,size_t>::iterator it2;
      for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
        it = row_QnumMdim.find(it2->first);
        _elemNum += it->second * it2->second;
      }
      return _elemNum;
    }

  template class Node<uni10_double64>;
  template class Node<uni10_complex128>;

}; /* namespace uni10 */
