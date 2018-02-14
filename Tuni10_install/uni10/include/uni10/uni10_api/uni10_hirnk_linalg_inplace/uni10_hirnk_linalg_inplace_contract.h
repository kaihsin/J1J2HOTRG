#ifndef __UNI10_HIGH_RANK_LINALG_INPLACE_CONTRACT_H__
#define __UNI10_HIGH_RANK_LINALG_INPLACE_CONTRACT_H__

#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{

  template<typename To, typename T, typename U>
    void contract( UniTensor<To>& Tc, UniTensor<T>& Ta, UniTensor<U>& Tb, bool fast, UNI10_INPLACE on);

  template<typename To, typename T, typename U>
    void contract( UniTensor<To>& Tc, UniTensor<T>& Ta, UniTensor<U>& Tb, bool fast, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );
      uni10_error_msg(!((*Ta.status) & (*Tb.status) & Ta.HAVEELEM), "%s" ,"Cannot perform contraction of two tensors before setting their elements.");

      if((void*)&Ta == (void*)&Tb || (void*)&Tc == (void*)&Tb){
        UniTensor<U> Ttmp(Tb);
        contract(Tc, Ta, Ttmp, fast, on);
        return ;
      }

      if((void*)&Tc == (void*)&Ta){
        UniTensor<T> Ttmp(Ta);
        contract(Tc, Ttmp, Tb, fast, on);
        return ;
      }

      if((*Ta.status) & Ta.HAVEBOND && (*Tb.status) & Ta.HAVEBOND){
        uni10_int AbondNum = Ta.bonds->size();
        uni10_int BbondNum = Tb.bonds->size();
        std::vector<uni10_int> oldLabelA = (*Ta.labels);
        std::vector<uni10_int> oldLabelB = (*Tb.labels);
        uni10_int oldRnumA = (*Ta.RBondNum);
        uni10_int oldRnumB = (*Tb.RBondNum);
        std::vector<uni10_int> newLabelA;
        std::vector<uni10_int> interLabel;
        std::vector<uni10_int> newLabelB;
        std::vector<uni10_int> markB(BbondNum, 0);
        std::vector<uni10_int> newLabelC;
        uni10_bool match;
        for(uni10_int a = 0; a < AbondNum; a++){
          match = false;
          for(uni10_int b = 0; b < BbondNum; b++)
            if((*Ta.labels)[a] == (*Tb.labels)[b]){
              markB[b] = 1;
              interLabel.push_back((*Ta.labels)[a]);
              newLabelB.push_back((*Tb.labels)[b]);

              uni10_error_msg(!( (*Ta.bonds)[a].dim() == (*Tb.bonds)[b].dim() ), "%s", 
                  "Cannot contract two bonds having different dimensions");

              match = true;
              break;
            }
          if(!match){
            newLabelA.push_back((*Ta.labels)[a]);
            newLabelC.push_back((*Ta.labels)[a]);
          }
        }
        for(uni10_uint64 a = 0; a < interLabel.size(); a++)
          newLabelA.push_back(interLabel[a]);
        for(uni10_int b = 0; b < BbondNum; b++)
          if(markB[b] == 0){
            newLabelB.push_back((*Tb.labels)[b]);
            newLabelC.push_back((*Tb.labels)[b]);
          }
        uni10_int conBond = interLabel.size();

        Ta = permute(Ta, newLabelA, AbondNum - conBond);
        Tb = permute(Tb, newLabelB, conBond);

        std::vector<Bond> cBonds;
        for(uni10_int i = 0; i < AbondNum - conBond; i++)
          cBonds.push_back((*Ta.bonds)[i]);
        for(uni10_int i = conBond; i < BbondNum; i++)
          cBonds.push_back((*Tb.bonds)[i]);

        Tc.assign(cBonds);
        if(cBonds.size())
          Tc.setLabel(newLabelC);
        Block<T> blockA;
        Block<U> blockB;
        Block<To> blockC;
        typename std::map<Qnum, Block<T> >::iterator it;
        typename std::map<Qnum, Block<U> >::iterator it2;
        for(it = Ta.blocks->begin() ; it != Ta.blocks->end(); it++){
          if((it2 = Tb.blocks->find(it->first)) != Tb.blocks->end()){
            blockA = it->second;
            blockB = it2->second;
            blockC = (*Tc.blocks)[it->first];
            uni10_error_msg(!(blockA.row() == blockC.row() && blockB.col() == blockC.col() && blockA.col() == blockB.row()), 
                "%s", "The dimensions the bonds to be contracted out are different.");
                       
            matrixDot(&blockA.elem_enforce(), &blockA.diag_enforce(), &blockB.elem_enforce(), &blockB.diag_enforce(), 
                &blockA.row_enforce(), &blockB.col_enforce(), &blockA.col_enforce(), &blockC.elem_enforce());
          }
        }
        (*Tc.status) |= Tc.HAVEELEM;

        if(conBond == 0){                     //Outer product

          uni10_int idx = 0;
          for(uni10_int i = 0; i < oldRnumA; i++){
            newLabelC[idx] = oldLabelA[i];
            idx++;
          }
          for(uni10_int i = 0; i < oldRnumB; i++){
            newLabelC[idx] = oldLabelB[i];
            idx++;
          }
          for(uni10_int i = oldRnumA; i < AbondNum; i++){
            newLabelC[idx] = oldLabelA[i];
            idx++;
          }
          for(uni10_int i = oldRnumB; i < BbondNum; i++){
            newLabelC[idx] = oldLabelB[i];
            idx++;
          }
          Tc = permute(Tc, newLabelC, oldRnumA + oldRnumB);
        }

        if(!fast){
          Ta = permute(Ta, oldLabelA, oldRnumA);
          Tb = permute(Tb, oldLabelB, oldRnumB);
        }
        return ;
      }
      else if((*Ta.status) & Ta.HAVEBOND)
        Tc = Ta * Tb[0];
      else if((*Tb.status) & Tb.HAVEBOND)
        Tc = Ta[0] * Tb;
      else
        Tc = UniTensor<To>(Ta[0] * Tb[0]);

    }

};

#endif
