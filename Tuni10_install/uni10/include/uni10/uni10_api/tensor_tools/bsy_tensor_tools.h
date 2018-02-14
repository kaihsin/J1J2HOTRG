#ifndef __UNI10_BSY_TENSOR_TOOLS_H__
#define __UNI10_BSY_TENSOR_TOOLS_H__

#include <stdio.h>

#include <set>

#include "uni10/uni10_api/UniTensor.h"

namespace uni10{

  namespace tensor_tools{

    //Function prototype.
    template<typename uni10_type>
      U_para<uni10_type>* init_para_bsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void copy_para_bsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para);

    template<typename uni10_type>
      void free_para_bsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void init_bsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      uni10_uint64 grouping_bsy(U_para<uni10_type>* _para);

    template <typename uni10_type>
      void initBlocks_bsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void setRawElem_bsy(U_para<uni10_type>* para, const uni10_type* rawElem);

    template <typename uni10_type>
      void putBlock_bsy(U_para<uni10_type>* para,const Qnum& qnum, const Block<uni10_type>& mat);

    template <typename uni10_type>
      void set_zero_bsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void randomize_bsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void transpose_bsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para);

    template <typename uni10_type>
      void dagger_bsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para);

    template <typename uni10_type>
      void conj_bsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para);

    template <typename uni10_type>
      void permute_bsy(const U_para<uni10_type>* T1_para, const std::vector<uni10_int32>& rsp_outin,
          U_para<uni10_type>* T2_para, uni10_bool inorder);

    template <typename uni10_type>
      void traceByRow_bsy(U_para<uni10_type>* Tout_para, const U_para<uni10_type>* Tin_para, uni10_int32 la, uni10_int32 lb );

    template <typename uni10_type>
      void addGate_bsy(U_para<uni10_type>* T1_para, const std::vector<_Swap>& swaps);

    template <typename uni10_type>
      uni10_type tensorAt_bsy(U_para<uni10_type>* T_para, const uni10_uint64* idxs);


    // Functions.
    template<typename uni10_type>
      U_para<uni10_type>* init_para_bsy(U_para<uni10_type>* para){

        U_para<uni10_type>* ptr;
        ptr      = new struct U_para<uni10_type>[1];
        ptr->bsy = new struct blk_sym_para<uni10_type>[1];
        ptr->check_status = 1;
        para = ptr;

        return para;

      }

    template<typename uni10_type>
      void copy_para_bsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para){

        uni10_error_msg(src_para->nsy != NULL, "%s", "Communication between nsy and bsy is developping!!");

        if(src_para!=NULL){
          *para->bsy = *src_para->bsy;
          typename std::map< Qnum, Block<uni10_type> >::const_iterator it2;
          typename std::map< const Block<uni10_type>* , Block<uni10_type>*> blkmap;
          typename std::map< Qnum, Block<uni10_type> >::iterator it = para->bsy->blocks.begin();
          for (; it != para->bsy->blocks.end(); it++ ){                       // blocks here is UniT.blocks
            it->second.elem_enforce().__elem = &(para->bsy->U_elem.__elem[it->second.elem_enforce().__elem - src_para->bsy->U_elem.__elem]);
            it2 = src_para->bsy->blocks.find(it->first);
            blkmap[&(it2->second)] = &(it->second);
          }

          if(src_para->bsy->status & UniTensor<uni10_type>::GET_HAVEBOND()){
            typename std::map<uni10_int32, Block<uni10_type>*>::iterator it = para->bsy->RQidx2Blk.begin();
            for(; it != para->bsy->RQidx2Blk.end(); it++)
              it->second = blkmap[it->second];
          }
        }

      }

    template<typename uni10_type>
      void free_para_bsy(U_para<uni10_type>* para){

        if(para != NULL){
          delete [] para->bsy;
          delete [] para;
        }

      }

    template<typename uni10_type>
      void init_bsy(U_para<uni10_type>* para){

        if(para->bsy->bonds.size()){
          para->bsy->U_elemNum = grouping_bsy(para);
          if(!(para->bsy->blocks.size() > 0)){      //No block in Tensor, Error!
            uni10_error_msg(true, "%s", "There is no symmetry block with the given bonds:\n");
            for(uni10_int32 b = 0; b < (uni10_int32)para->bsy->bonds.size(); b++)
              std::cout<<"    "<<para->bsy->bonds[b];
          }

          para->bsy->labels.assign(para->bsy->bonds.size(), 0);
          for(uni10_int32 b = 0; b < (uni10_int32)para->bsy->bonds.size(); b++)
            para->bsy->labels[b] = b;
          //para->bsy->status |= UniTensor<uni10_type>::GET_HAVEBOND();

        }
        else{
          Qnum q0(0);
          para->bsy->blocks[q0] = Block<uni10_type>(1, 1);
          para->bsy->RBondNum = 0;
          para->bsy->RQdim = 0;
          para->bsy->CQdim = 0;
          para->bsy->U_elemNum = 1;
          //para->bsy->status |= UniTensor<uni10_type>::GET_HAVEELEM();
        }

        para->bsy->U_elem.init(1, para->bsy->U_elemNum, false);
        initBlocks_bsy(para);

      }

    template<typename uni10_type>
      uni10_uint64 grouping_bsy(U_para<uni10_type>* para){

        para->bsy->blocks.clear();
        uni10_int32 row_bondNum = 0;
        uni10_int32 col_bondNum = 0;
        para->bsy->RQdim = 1;
        para->bsy->CQdim = 1;
        uni10_bool IN_BONDS_BEFORE_OUT_BONDS = true;

        for(uni10_uint64 i = 0; i < para->bsy->bonds.size(); i++){

          if(para->bsy->bonds[i].type() == BD_IN){
            uni10_error_msg(!(IN_BONDS_BEFORE_OUT_BONDS == true), 
                "%s","Error in the input bond array: BD_OUT bonds must be placed after all BD_IN bonds.");

            para->bsy->RQdim *= para->bsy->bonds[i].const_getQnums().size();
            row_bondNum++;
          }
          else{
            para->bsy->CQdim *=  para->bsy->bonds[i].const_getQnums().size();
            col_bondNum++;
            IN_BONDS_BEFORE_OUT_BONDS = false;
          }
        }
        para->bsy->RBondNum = row_bondNum;

        std::map<Qnum,uni10_uint64> row_QnumMdim;
        std::vector<uni10_int32> row_offs(row_bondNum, 0);
        std::map<Qnum,std::vector<uni10_int32> > row_Qnum2Qidx;
        Qnum qnum;

        uni10_uint64 dim;
        uni10_int32 boff = 0;

        std::vector<uni10_uint64>tmpRQidx2Dim(para->bsy->RQdim, 1);
        std::vector<uni10_uint64>tmpCQidx2Dim(para->bsy->CQdim, 1);
        std::vector<uni10_uint64>tmpRQidx2Off(para->bsy->RQdim, 0);
        std::vector<uni10_uint64>tmpCQidx2Off(para->bsy->CQdim, 0);

        if(row_bondNum){
          while(1){
            qnum.assign();
            dim = 1;
            for(uni10_int32 b = 0; b < row_bondNum; b++){
              qnum = qnum * para->bsy->bonds[b].const_getQnums()[row_offs[b]];
              dim *= para->bsy->bonds[b].const_getQdegs()[row_offs[b]];
            }
            if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
              tmpRQidx2Off[boff] = row_QnumMdim[qnum];
              tmpRQidx2Dim[boff] = dim;
              row_QnumMdim[qnum] += dim;
            }
            else{
              tmpRQidx2Off[boff] = 0;
              tmpRQidx2Dim[boff] = dim;
              row_QnumMdim[qnum] = dim;
            }
            row_Qnum2Qidx[qnum].push_back(boff);
            boff++;
            uni10_int32 bidx;
            for(bidx = row_bondNum - 1; bidx >= 0; bidx--){
              row_offs[bidx]++;
              if(row_offs[bidx] < para->bsy->bonds[bidx].const_getQnums().size())
                break;
              else
                row_offs[bidx] = 0;
            }
            if(bidx < 0)  //run over all row_bond offsets
              break;
          }
        }
        else{
          qnum.assign();
          row_QnumMdim[qnum] = 1;
          row_Qnum2Qidx[qnum].push_back(0);
        }
        std::map<Qnum,uni10_uint64> col_QnumMdim;
        std::vector<uni10_int32> col_offs(col_bondNum, 0);
        std::map<Qnum,std::vector<uni10_int32> > col_Qnum2Qidx;
        boff = 0;
        if(col_bondNum){
          while(1){
            qnum.assign();
            dim = 1;
            for(uni10_int32 b = 0; b < col_bondNum; b++){
              qnum = qnum * para->bsy->bonds[b + row_bondNum].const_getQnums()[col_offs[b]];
              dim *= para->bsy->bonds[b + row_bondNum].const_getQdegs()[col_offs[b]];
            }
            if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
              if(col_QnumMdim.find(qnum) != col_QnumMdim.end()){
                tmpCQidx2Off[boff] = col_QnumMdim[qnum];
                tmpCQidx2Dim[boff] = dim;
                col_QnumMdim[qnum] += dim;
              }
              else{
                tmpCQidx2Off[boff] = 0;
                tmpCQidx2Dim[boff] = dim;
                col_QnumMdim[qnum] = dim;
              }
              col_Qnum2Qidx[qnum].push_back(boff);
            }
            boff++;
            uni10_int32 bidx;
            for(bidx = col_bondNum - 1; bidx >= 0; bidx--){
              col_offs[bidx]++;
              if(col_offs[bidx] < para->bsy->bonds[bidx + row_bondNum].const_getQnums().size())
                break;
              else
                col_offs[bidx] = 0;
            }
            if(bidx < 0)  //run over all row_bond offsets
              break;
          }
        }
        else{
          qnum.assign();
          if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
            col_QnumMdim[qnum] = 1;
            col_Qnum2Qidx[qnum].push_back(0);
          }
        }

        std::map<Qnum,uni10_uint64>::iterator it;
        std::map<Qnum,uni10_uint64>::iterator it2;
        std::set<uni10_int32> Qidx;
        uni10_int32 qidx;
        uni10_uint64 off = 0;
        for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
          it = row_QnumMdim.find(it2->first);
          Block<uni10_type> blk(it->second, it2->second); // blk(Rnum, Cnum);
          off += blk.row() * blk.col();
          para->bsy->blocks[it->first] = blk;
          Block<uni10_type>* blkptr = &(para->bsy->blocks[it->first]);
          std::vector<uni10_int32>& tmpRQidx = row_Qnum2Qidx[it->first];
          std::vector<uni10_int32>& tmpCQidx = col_Qnum2Qidx[it->first];
          for(uni10_uint64 i = 0; i < tmpRQidx.size(); i++){
            para->bsy->RQidx2Blk[tmpRQidx[i]] = blkptr;
            for(uni10_uint64 j = 0; j < tmpCQidx.size(); j++){
              para->bsy->RQidx2Dim[tmpRQidx[i]] = tmpRQidx2Dim[tmpRQidx[i]];
              para->bsy->RQidx2Off[tmpRQidx[i]] = tmpRQidx2Off[tmpRQidx[i]];
              para->bsy->CQidx2Dim[tmpCQidx[j]] = tmpCQidx2Dim[tmpCQidx[j]];
              para->bsy->CQidx2Off[tmpCQidx[j]] = tmpCQidx2Off[tmpCQidx[j]];
              qidx = tmpRQidx[i] * para->bsy->CQdim + tmpCQidx[j];
              Qidx.insert(qidx);
            }
          }
        }
        uni10_uint64 elemEnc = 0;
        for(std::map<uni10_int32, uni10_uint64>::iterator itr = para->bsy->RQidx2Dim.begin(); itr != para->bsy->RQidx2Dim.end(); itr++)
          for(std::map<uni10_int32, uni10_uint64>::iterator itc = para->bsy->CQidx2Dim.begin(); itc != para->bsy->CQidx2Dim.end(); itc++){
            qidx = itr->first * para->bsy->CQdim + itc->first;
            if(Qidx.find(qidx) != Qidx.end()){
              para->bsy->QidxEnc[qidx] = elemEnc;
              elemEnc += para->bsy->RQidx2Dim[itr->first] * para->bsy->CQidx2Dim[itc->first];
            }
          }
        return off;

      }

    template <typename uni10_type>
      void initBlocks_bsy(U_para<uni10_type>* para){

        uni10_uint64 offset = 0;
        typename std::map< Qnum, Block<uni10_type> >::iterator it = para->bsy->blocks.begin();
        for(; it != para->bsy->blocks.end(); it++ ){
          it->second.elem_enforce().__elem = &(para->bsy->U_elem.__elem[offset]);
          offset += it->second.row_enforce() * it->second.col_enforce();
        }

      }

    template <typename uni10_type>
      void putBlock_bsy(U_para<uni10_type>* para,const Qnum& qnum, const Block<uni10_type>& mat){

        typename std::map<Qnum, Block<uni10_type> >::iterator it;

        if(!((it = para->bsy->blocks.find(qnum)) != para->bsy->blocks.end())){
          uni10_error_msg(true, "%s", "There is no block with the given quantum number ");
          std::cout<<qnum;
        }

        uni10_error_msg(!(mat.row() == it->second.row_enforce() && mat.col() == it->second.col_enforce()), "%s", 
            "The dimension of input matrix does not match for the dimension of the block with quantum number \n  Hint: Use Matrix::resize(int, int)");

        if(mat.getElem() != it->second.getElem()){
          if(mat.isDiag()){
            setDiag(&it->second.elem_enforce(), &mat.const_elem_enforce(), &it->second.row_enforce(), &it->second.col_enforce());
          }
          else
            it->second.elem_enforce().copy(0, mat.const_elem_enforce(), it->second.row_enforce() * it->second.col_enforce() );

        }

        para->bsy->status |= UniTensor<uni10_type>::GET_HAVEELEM();

      }

    template <typename uni10_type>
      void set_zero_bsy(U_para<uni10_type>* para){

        para->bsy->U_elem.set_zeros();

      }

    template <typename uni10_type>
      void randomize_bsy(U_para<uni10_type>* para){

        para = NULL;
        uni10_error_msg(true, "%s", "Developping");

      }

    template <typename uni10_type>
      void setRawElem_bsy(U_para<uni10_type>* para, const uni10_type* rawElem){

        uni10_error_msg((para->bsy->status & UniTensor<uni10_type>::GET_HAVEBOND()) == 0, 
            "%s", "Setting elements to a tensor without bonds is not supported." );

        uni10_int32 bondNum = para->bsy->bonds.size();
        std::vector<uni10_int32> Q_idxs(bondNum, 0);
        std::vector<uni10_int32> Q_Bdims(bondNum, 0);
        std::vector<uni10_int32> sB_idxs(bondNum, 0);
        std::vector<uni10_int32> sB_sBdims(bondNum, 0);
        std::vector<uni10_int32> rAcc(bondNum, 1);
        for(uni10_int32 b = 0; b < bondNum; b++)
          Q_Bdims[b] = para->bsy->bonds[b].const_getQnums().size();
        for(uni10_int32 b = bondNum - 1; b > 0; b--)
          rAcc[b - 1] = rAcc[b] * para->bsy->bonds[b].dim();
        uni10_int32 Q_off;
        uni10_int32 tmp;
        uni10_int32 RQoff, CQoff;
        uni10_uint64 sB_r, sB_c;        //sub-block of a Qidx
        uni10_uint64 sB_rDim, sB_cDim;  //sub-block of a Qidx
        uni10_uint64 B_cDim;
        uni10_uint64 E_off;
        uni10_int32 R_off;
        uni10_type* work = para->bsy->U_elem.__elem;

        for(std::map<uni10_int32, uni10_uint64>::iterator it = para->bsy->QidxEnc.begin(); it != para->bsy->QidxEnc.end(); it++){
          Q_off = it->first;
          tmp = Q_off;
          for(uni10_int32 b = bondNum - 1; b >= 0; b--){
            Q_idxs[b] = tmp % Q_Bdims[b];
            tmp /= Q_Bdims[b];
          }
          R_off = 0;
          for(uni10_int32 b = 0; b < bondNum; b++){
            R_off += rAcc[b] * para->bsy->bonds[b].const_getOffsets()[Q_idxs[b]];
            sB_sBdims[b] = para->bsy->bonds[b].const_getQdegs()[Q_idxs[b]];
          }
          RQoff = Q_off / para->bsy->CQdim;
          CQoff = Q_off % para->bsy->CQdim;
          B_cDim = para->bsy->RQidx2Blk[RQoff]->col();

          E_off = (para->bsy->RQidx2Blk[RQoff]->elem_enforce().__elem - para->bsy->U_elem.__elem) + (para->bsy->RQidx2Off[RQoff] * B_cDim) + para->bsy->CQidx2Off[CQoff];

          sB_rDim = para->bsy->RQidx2Dim[RQoff];
          sB_cDim = para->bsy->CQidx2Dim[CQoff];
          sB_idxs.assign(bondNum, 0);
          for(sB_r = 0; sB_r < sB_rDim; sB_r++)
            for(sB_c = 0; sB_c < sB_cDim; sB_c++){
              work[E_off + (sB_r * B_cDim) + sB_c] = rawElem[R_off];
              for(uni10_int32 bend = bondNum - 1; bend >= 0; bend--){
                sB_idxs[bend]++;
                if(sB_idxs[bend] < sB_sBdims[bend]){
                  R_off += rAcc[bend];
                  break;
                }
                else{
                  R_off -= rAcc[bend] * (sB_idxs[bend] - 1);
                  sB_idxs[bend] = 0;
                }
              }
            }
        }

        para->bsy->status |= UniTensor<uni10_type>::GET_HAVEELEM();

      }

    template <typename uni10_type>
      void transpose_bsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para){

        typename std::map<Qnum, Block<uni10_type> >::iterator it_in;
        typename std::map<Qnum, Block<uni10_type> >::iterator it_out;
        UELEM(uni10_elem, _package, _type)<uni10_type>* elem_in;
        UELEM(uni10_elem, _package, _type)<uni10_type>* elem_out;
        uni10_uint64 Rnum, Cnum;
        for ( it_in = src_para->bsy->blocks.begin() ; it_in != src_para->bsy->blocks.end(); it_in++ ){
          it_out = para->bsy->blocks.find((it_in->first));
          Rnum = it_in->second.row_enforce();
          Cnum = it_in->second.col_enforce();
          elem_in  = &(it_in->second.elem_enforce());
          elem_out = &(it_out->second.elem_enforce());
          setTranspose(elem_in, &Rnum, &Cnum, elem_out);
        }

      }


    template <typename uni10_type>
      void dagger_bsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para){

        typename std::map<Qnum, Block<uni10_type> >::iterator it_in;
        typename std::map<Qnum, Block<uni10_type> >::iterator it_out;
        UELEM(uni10_elem, _package, _type)<uni10_type>* elem_in;
        UELEM(uni10_elem, _package, _type)<uni10_type>* elem_out;
        uni10_uint64 Rnum, Cnum;
        for ( it_in = src_para->bsy->blocks.begin() ; it_in != src_para->bsy->blocks.end(); it_in++ ){
          it_out = para->bsy->blocks.find((it_in->first));
          Rnum = it_in->second.row_enforce();
          Cnum = it_in->second.col_enforce();
          elem_in  = &(it_in->second.elem_enforce());
          elem_out = &(it_out->second.elem_enforce());
          setDagger(elem_in, &Rnum, &Cnum, elem_out);
        }

      }

    template <typename uni10_type>
      void conj_bsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para){

        setConjugate(&src_para->bsy->U_elem, &src_para->bsy->U_elemNum, &para->bsy->U_elem);

      }

    template <typename uni10_type>
      void permute_bsy(const U_para<uni10_type>* T1_para, const std::vector<uni10_int32>& rsp_outin,
          U_para<uni10_type>* T2_para, uni10_bool inorder){

        uni10_uint64 bondNum = T1_para->bsy->bonds.size();
        uni10_double64 sign = 1.0;
        //For Fermionic system
        std::vector<_Swap> swaps;
        if(Qnum::isFermionic()){
          std::vector<uni10_int32> inLabelF(bondNum);
          std::vector<uni10_int32> outLabelF(bondNum);
          std::vector<uni10_int32> ordF(bondNum);

          for(uni10_int32 b = 0; b < T1_para->bsy->RBondNum; b++){
            inLabelF[b] = T1_para->bsy->labels[b];
            ordF[b] = b;
          }
          for(uni10_int32 b = 0; b < T2_para->bsy->RBondNum; b++)
            outLabelF[b] = T2_para->bsy->labels[b];
          for(uni10_int32 b = bondNum - 1; b >= T1_para->bsy->RBondNum; b--){
            ordF[b] = bondNum - b + T1_para->bsy->RBondNum - 1;
            inLabelF[ordF[b]] = T1_para->bsy->labels[b];
          }
          for(uni10_int32 b = bondNum - 1; b >= T2_para->bsy->RBondNum; b--)
            outLabelF[bondNum - b + T2_para->bsy->RBondNum - 1] = T2_para->bsy->labels[b];

          std::vector<uni10_int32> rspF_outin(bondNum);
          for(uni10_uint64 i = 0; i < bondNum; i++)
            for(uni10_uint64 j = 0; j < bondNum; j++)
              if(inLabelF[i] == outLabelF[j])
                rspF_outin[j] = i;
          swaps = recSwap(rspF_outin, ordF);
        }

        //End Fermionic system
        std::vector<uni10_int32> Qin_idxs(bondNum, 0);
        std::vector<uni10_int32> Qot_idxs(bondNum, 0);
        uni10_int32 Qin_off, Qot_off;
        uni10_int32 tmp;
        uni10_int32 Qin_RQoff, Qin_CQoff;
        uni10_int32 Qot_CQoff, Qot_RQoff;
        uni10_uint64 sBin_r, sBin_c;	//sub-block of a Qidx
        uni10_uint64 sBin_rDim, sBin_cDim;	//sub-block of a Qidx
        uni10_uint64 sBot_cDim;	//sub-block of a Qidx
        uni10_uint64 sBot_r, sBot_c;
        uni10_uint64 Bin_cDim, Bot_cDim;
        uni10_type* Ein_ptr;
        uni10_type* Eot_ptr;
        std::vector<uni10_int32> sBin_idxs(bondNum, 0);
        std::vector<uni10_int32> sBin_sBdims(bondNum, 0);
        std::vector<uni10_int32> Qot_acc(bondNum, 1);
        std::vector<uni10_int32> sBot_acc(bondNum, 1);
        for(uni10_int32 b = bondNum	- 1; b > 0; b--)
          Qot_acc[b - 1] = Qot_acc[b] * T2_para->bsy->bonds[b].const_getQnums().size();

        for(std::map<uni10_int32, uni10_uint64>::iterator it = T1_para->bsy->QidxEnc.begin(); it != T1_para->bsy->QidxEnc.end(); it++){
          Qin_off = it->first;
          tmp = Qin_off;
          uni10_int32 qdim;
          for(uni10_int32 b = bondNum - 1; b >= 0; b--){
            qdim = T1_para->bsy->bonds[b].const_getQnums().size();
            Qin_idxs[b] = tmp % qdim;
            sBin_sBdims[b] = T1_para->bsy->bonds[b].const_getQdegs()[Qin_idxs[b]];
            tmp /= qdim;
          }
          Qot_off = 0;
          for(uni10_uint64 b = 0; b < bondNum; b++){
            Qot_idxs[b] = Qin_idxs[rsp_outin[b]];
            Qot_off += Qot_idxs[b] * Qot_acc[b];
          }
          for(uni10_uint64 b = bondNum - 1; b > 0; b--)
            sBot_acc[rsp_outin[b-1]] = sBot_acc[rsp_outin[b]] * T1_para->bsy->bonds[rsp_outin[b]].const_getQdegs()[Qot_idxs[b]];

          Qin_RQoff = Qin_off / T1_para->bsy->CQdim;
          Qin_CQoff = Qin_off % T1_para->bsy->CQdim;
          Qot_RQoff = Qot_off / T2_para->bsy->CQdim;
          Qot_CQoff = Qot_off % T2_para->bsy->CQdim;
          Bin_cDim = T1_para->bsy->RQidx2Blk[Qin_RQoff]->col();
          Bot_cDim = T2_para->bsy->RQidx2Blk[Qot_RQoff]->col();
          Ein_ptr = T1_para->bsy->RQidx2Blk[Qin_RQoff]->getElem() + (T1_para->bsy->RQidx2Off[Qin_RQoff] * Bin_cDim) + T1_para->bsy->CQidx2Off[Qin_CQoff];
          Eot_ptr = T2_para->bsy->RQidx2Blk[Qot_RQoff]->getElem() + (T2_para->bsy->RQidx2Off[Qot_RQoff] * Bot_cDim) + T2_para->bsy->CQidx2Off[Qot_CQoff];
          sBin_rDim = T1_para->bsy->RQidx2Dim[Qin_RQoff];
          sBin_cDim = T1_para->bsy->CQidx2Dim[Qin_CQoff];
          sBot_cDim = T2_para->bsy->CQidx2Dim[Qot_CQoff];
          uni10_int32 cnt_ot = 0;
          sBin_idxs.assign(bondNum, 0);
          if(Qnum::isFermionic()){
            uni10_int32 sign01 = 0;
            for(uni10_uint64 i = 0; i < swaps.size(); i++)
              sign01 ^= (T1_para->bsy->bonds[swaps[i].b1].const_getQnums()[Qin_idxs[swaps[i].b1]].prtF() & T1_para->bsy->bonds[swaps[i].b2].const_getQnums()[Qin_idxs[swaps[i].b2]].prtF());
            sign = sign01 ? -1.0 : 1.0;
          }
          for(sBin_r = 0; sBin_r < sBin_rDim; sBin_r++)
            for(sBin_c = 0; sBin_c < sBin_cDim; sBin_c++){
              sBot_r = cnt_ot / sBot_cDim;
              sBot_c = cnt_ot % sBot_cDim;
              Eot_ptr[(sBot_r * Bot_cDim) + sBot_c] = sign * Ein_ptr[(sBin_r * Bin_cDim) + sBin_c];
              for(uni10_int32 bend = bondNum - 1; bend >= 0; bend--){
                sBin_idxs[bend]++;
                if(sBin_idxs[bend] < sBin_sBdims[bend]){
                  cnt_ot += sBot_acc[bend];
                  break;
                }
                else{
                  cnt_ot -= sBot_acc[bend] * (sBin_idxs[bend] - 1);
                  sBin_idxs[bend] = 0;
                }
              }
            }
        }

      }

    template <typename uni10_type>
      void addGate_bsy(U_para<uni10_type>* T1_para, const std::vector<_Swap>& swaps){

        uni10_error_msg(true, "%s", "Developping");
        T1_para = NULL;
        int a = swaps.size();
        a = 0;

      }

    template <typename uni10_type>
      void traceByRow_bsy(U_para<uni10_type>* Tout_para, const U_para<uni10_type>* Tin_para, uni10_int32 ia, uni10_int32 ib ){

        uni10_int32 bondNum = Tin_para->bsy->bonds.size();
        std::vector<uni10_int32> Q_acc(bondNum, 1);
        for(uni10_int32 b = bondNum - 1; b > 0; b--)
          Q_acc[b - 1] = Q_acc[b] * Tin_para->bsy->bonds[b].const_getQnums().size();

        uni10_int32 tQdim = Tin_para->bsy->bonds[ia].const_getQnums().size();
        uni10_error_msg(  !(tQdim == Tin_para->bsy->bonds[ib].const_getQnums().size()), "%s",  "The bonds of the given two labels does not match for trace.");

        Qnum q0(0, PRT_EVEN);
        for(uni10_int32 q = 0; q < tQdim; q++){
          uni10_error_msg(!((Tin_para->bsy->bonds[ia].const_getQnums()[q] * Tin_para->bsy->bonds[ib].const_getQnums()[q] == q0) 
                && (Tin_para->bsy->bonds[ia].const_getQdegs()[q] == Tin_para->bsy->bonds[ib].const_getQdegs()[q]))
              , "%s", "The bonds of the given two labels does not match for trace.");
        }

        uni10_int32 tBnum = Tout_para->bsy->bonds.size();
        std::vector<uni10_int32> Qt_Bdims(tBnum, 0);
        for(uni10_int32 b = 0; b < tBnum; b++)
          Qt_Bdims[b] = Tout_para->bsy->bonds[b].const_getQnums().size();

        uni10_int32 Qt_off;
        uni10_int32 Q_off;
        uni10_int32 Qt_RQoff, Qt_CQoff;
        uni10_int32 Q_RQoff, Q_CQoff;
        uni10_uint64 sBt_rDim, sBt_cDim;  //sub-block of a Qidx of Tt
        uni10_uint64 Bt_cDim;
        uni10_type* Et_ptr;

        std::vector<uni10_type*> E_offs(tQdim);
        std::vector<uni10_uint64> B_cDims(tQdim);
        uni10_int32 tQdim2 = tQdim * tQdim;
        uni10_int32 Qenc = Q_acc[ia] + Q_acc[ib];

        typename std::map<uni10_int32, uni10_uint64>::iterator it = Tout_para->bsy->QidxEnc.begin();

        for(; it != Tout_para->bsy->QidxEnc.end(); it++){

          Qt_off = it->first;
          Qt_RQoff = Qt_off / Tout_para->bsy->CQdim;
          Qt_CQoff = Qt_off % Tout_para->bsy->CQdim;
          Bt_cDim = Tout_para->bsy->RQidx2Blk[Qt_RQoff]->col_enforce();
          Et_ptr = Tout_para->bsy->RQidx2Blk[Qt_RQoff]->elem_enforce().__elem + (Tout_para->bsy->RQidx2Off[Qt_RQoff] * Bt_cDim) + Tout_para->bsy->CQidx2Off[Qt_CQoff];
          sBt_rDim = Tout_para->bsy->RQidx2Dim[Qt_RQoff];
          sBt_cDim = Tout_para->bsy->CQidx2Dim[Qt_CQoff];

          for(uni10_int32 q = 0; q < tQdim; q++){
            Q_off = Qt_off * tQdim2 + q * Qenc;
            Q_RQoff = Q_off / Tin_para->bsy->CQdim;
            Q_CQoff = Q_off % Tin_para->bsy->CQdim;
            B_cDims[q] = Tin_para->bsy->RQidx2Blk[Q_RQoff]->col_enforce();
            E_offs[q] = Tin_para->bsy->RQidx2Blk[Q_RQoff]->elem_enforce().__elem + (Tin_para->bsy->RQidx2Off[Q_RQoff] * B_cDims[q]) + Tin_para->bsy->CQidx2Off[Q_CQoff];
          }
          uni10_int32 tQdeg, sB_c_off;
          uni10_type trVal;
          for(uni10_uint64 sB_r = 0; sB_r < sBt_rDim; sB_r++)
            for(uni10_uint64 sB_c = 0; sB_c < sBt_cDim; sB_c++){
              trVal = 0;
              for(uni10_int32 q = 0; q < tQdim; q++){
                tQdeg = Tin_para->bsy->bonds[ia].const_getQdegs()[q];
                sB_c_off = sB_c * (tQdeg * tQdeg);
                for(uni10_int32 t = 0; t < tQdeg; t++){
                  trVal += E_offs[q][(sB_r * B_cDims[q]) + sB_c_off + t * (tQdeg + 1)];
                }
              }
              Et_ptr[sB_r * Bt_cDim + sB_c] = trVal;
            }
        }

      }

    template <typename uni10_type>
      uni10_type tensorAt_bsy(U_para<uni10_type>* T_para, const uni10_uint64* idxs){

        uni10_int32 bondNum = T_para->bsy->bonds.size();
        std::vector<uni10_int32> Qidxs(bondNum, 0);
        for(uni10_int32 b = 0; b < bondNum; b++){
          uni10_error_msg(!(idxs[b] < (uni10_uint64)T_para->bsy->bonds[b].dim()), "%s", 
              "The input indices are out of range.");
          for(uni10_int32 q = T_para->bsy->bonds[b].const_getOffsets().size() - 1; q >= 0; q--){
            if(idxs[b] < (T_para->bsy->bonds)[b].const_getOffsets()[q])
              continue;
            Qidxs[b] = q;
            break;
          }
        }

        std::vector<uni10_int32> Q_acc(bondNum, 1);
        for(uni10_int32 b = bondNum	- 1; b > 0; b--)
          Q_acc[b - 1] = Q_acc[b] * (T_para->bsy->bonds)[b].const_getQnums().size();
        uni10_int32 Qoff = 0;
        for(uni10_int32 b = 0; b < bondNum; b++)
          Qoff += Q_acc[b] * Qidxs[b];

        if(T_para->bsy->QidxEnc.find(Qoff) != T_para->bsy->QidxEnc.end()){
          uni10_int32 Q_RQoff = Qoff / T_para->bsy->CQdim;
          uni10_int32 Q_CQoff = Qoff % T_para->bsy->CQdim;
          Block<uni10_type>* blk = T_para->bsy->RQidx2Blk.find(Q_RQoff)->second;
          uni10_uint64 B_cDim = blk->col_enforce();
          uni10_uint64 sB_cDim = T_para->bsy->CQidx2Dim.find(Q_CQoff)->second;
          uni10_uint64 blkRoff = T_para->bsy->RQidx2Off.find(Q_RQoff)->second;
          uni10_uint64 blkCoff = T_para->bsy->CQidx2Off.find(Q_CQoff)->second;

          uni10_type* boff = &blk->elem_enforce()[0] + (blkRoff * B_cDim) + blkCoff;

          uni10_int32 cnt = 0;
          std::vector<uni10_int32> D_acc(bondNum, 1);
          for(uni10_int32 b = bondNum	- 1; b > 0; b--)
            D_acc[b - 1] = D_acc[b] * T_para->bsy->bonds[b].const_getQdegs()[Qidxs[b]];
          for(uni10_int32 b = 0; b < bondNum; b++)
            cnt += (idxs[b] - T_para->bsy->bonds[b].const_getOffsets()[Qidxs[b]]) * D_acc[b];
          return boff[(cnt / sB_cDim) * B_cDim + cnt % sB_cDim];
        }
        else{
          return 0.0;
        }

      }

  };

};

#endif
