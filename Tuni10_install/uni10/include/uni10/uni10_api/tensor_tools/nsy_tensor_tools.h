#ifndef __UNI10_NSY_TENSOR_TOOLS_H__
#define __UNI10_NSY_TENSOR_TOOLS_H__

#include <stdio.h>

#include "uni10/uni10_api/rng.h"
#include "uni10/uni10_api/UniTensor.h"

namespace uni10{
  
  namespace tensor_tools{

    //Function prototype.
    template<typename uni10_type>
      U_para<uni10_type>* init_para_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void copy_para_nsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para);

    template<typename uni10_type>
      void free_para_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      void init_nsy(U_para<uni10_type>* para);

    template<typename uni10_type>
      uni10_uint64 grouping_nsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void initBlocks_nsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void setRawElem_nsy(U_para<uni10_type>* para, const uni10_type* rawElem);

    template <typename uni10_type>
      void putBlock_nsy(U_para<uni10_type>* para,const Qnum& qnum, const Block<uni10_type>& mat);

    template <typename uni10_type>
      void set_zero_nsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void randomize_nsy(U_para<uni10_type>* para);

    template <typename uni10_type>
      void transpose_nsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para);

    template <typename uni10_type>
      void dagger_nsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para);

    template <typename uni10_type>
      void conj_nsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para);

    template <typename uni10_type>
      void permute_nsy(const U_para<uni10_type>* T1_para, const std::vector<uni10_int32>& rsp_outin,
        U_para<uni10_type>* T2_para, uni10_bool inorder);

    template <typename uni10_type>
      void traceByRow_nsy(U_para<uni10_type>* Tout_para, const U_para<uni10_type>* Tin_para, uni10_int32 la, uni10_int32 lb );

    template <typename uni10_type>
      void addGate_nsy(U_para<uni10_type>* T1_para, const std::vector<_Swap>& swaps);

    template <typename uni10_type>
      uni10_type tensorAt_nsy(U_para<uni10_type>* T_para, const uni10_uint64* idxs);


    //Functions.
    template<typename uni10_type>
       U_para<uni10_type>* init_para_nsy(U_para<uni10_type>* para){

        U_para<uni10_type>* ptr;
        ptr      = new struct U_para<uni10_type>[1];
        ptr->nsy = new struct no_sym_para<uni10_type>[1];
        ptr->check_status = 1;
        para = ptr;

        return para;
        
      }

    template<typename uni10_type>
      void copy_para_nsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para){

        if(src_para != NULL)
          *para->nsy = *src_para->nsy;

      }

    template<typename uni10_type>
       void free_para_nsy(U_para<uni10_type>* para){
        
         if(para!=NULL){
           delete [] para->nsy;
           delete [] para;
         }
        
      }

    template<typename uni10_type>
      void init_nsy(U_para<uni10_type>* para){

        if(para->nsy->bonds.size()){
          para->nsy->U_elemNum = grouping_nsy(para);
          if(!(para->nsy->blocks.size() > 0)){ //No block in Tensor, Error!
            uni10_error_msg(true, "%s", "There is no symmetry block with the given bonds:\n");
            for(uni10_int32 b = 0; b < (uni10_int32)para->nsy->bonds.size(); b++)
              std::cout<<"    "<<para->nsy->bonds[b];
          }

          para->nsy->labels.assign(para->nsy->bonds.size(), 0);
          for(uni10_int32 b = 0; b < (uni10_int32)para->nsy->bonds.size(); b++)
            para->nsy->labels[b] = b;
          //para->nsy->status |= UniTensor<uni10_type>::GET_HAVEBOND();

        }
        else{
          Qnum q0(0);
          para->nsy->blocks[q0] = Block<uni10_type>(1, 1);
          para->nsy->RBondNum = 0;
          para->nsy->RQdim = 0;
          para->nsy->CQdim = 0;
          para->nsy->U_elemNum = 1;
          //para->nsy->status |= UniTensor<uni10_type>::GET_HAVEELEM();
        }

        para->nsy->U_elem.init(1, para->nsy->U_elemNum, false);
        //para->nsy->U_elem.print_elem(para->nsy->RQdim, para->nsy->U_elemNum, false)
        initBlocks_nsy(para);

      }

    template<typename uni10_type>
      uni10_uint64 grouping_nsy(U_para<uni10_type>* para){

        para->nsy->blocks.clear();
        Qnum q0(0);
        uni10_int32   row_bondNum = 0;
        uni10_int32   col_bondNum = 0;
        para->nsy->RQdim = 1;
        para->nsy->CQdim = 1;
        uni10_bool IN_BONDS_BEFORE_OUT_BONDS = true;
        for(uni10_int32 i = 0; i < (uni10_int32)para->nsy->bonds.size(); i++){
          uni10_error_msg( para->nsy->bonds[i].const_getQdegs().size() > 1 || para->nsy->bonds[i].const_getQnums().size() > 1
              , "%s", "This UniTensor has symmetry!!");
          if(para->nsy->bonds[i].type() == BD_IN){
            uni10_error_msg(!(IN_BONDS_BEFORE_OUT_BONDS == true), 
                "%s","Error in the input bond array: BD_OUT bonds must be placed after all BD_IN bonds.");
            para->nsy->RQdim *= para->nsy->bonds[i].const_getQdegs()[0];
            row_bondNum++;
          }
          else{
            para->nsy->CQdim *= para->nsy->bonds[i].const_getQdegs()[0];
            col_bondNum++;
            IN_BONDS_BEFORE_OUT_BONDS = false;
          }
        }
        para->nsy->RBondNum = row_bondNum;
        para->nsy->blocks[q0] = Block<uni10_type>(para->nsy->RQdim, para->nsy->CQdim);

        return para->nsy->RQdim*para->nsy->CQdim;

      }

    template <typename uni10_type>
      void initBlocks_nsy(U_para<uni10_type>* para){
        uni10_uint64 offset = 0;
        typename std::map< Qnum, Block<uni10_type> >::iterator it = para->nsy->blocks.begin();
        for(; it != para->nsy->blocks.end(); it++ ){
          it->second.elem_enforce().__elem = &(para->nsy->U_elem.__elem[offset]);
          offset += it->second.row_enforce() * it->second.col_enforce();
        }
      }

    template <typename uni10_type>
      void setRawElem_nsy(U_para<uni10_type>* para, const uni10_type* rawElem){
      
        para->nsy->U_elem.setElem(rawElem);

      }

    template <typename uni10_type>
      void putBlock_nsy(U_para<uni10_type>* para,const Qnum& qnum, const Block<uni10_type>& mat){

        typename std::map<Qnum, Block<uni10_type> >::iterator it;

        if(!((it = para->nsy->blocks.find(qnum)) != para->nsy->blocks.end())){
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

      }

    template <typename uni10_type>
      void set_zero_nsy(U_para<uni10_type>* para){

        para->nsy->U_elem.set_zeros();

      }

    template <typename uni10_type>
      void randomize_nsy(U_para<uni10_type>* para){

        uni10_error_msg(true, "%s", "Developping!!!\n");

      }

    template <typename uni10_type>
      void transpose_nsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para){

        typename std::map<Qnum, Block<uni10_type> >::iterator it_in = src_para->nsy->blocks.begin();
        typename std::map<Qnum, Block<uni10_type> >::iterator it_out = para->nsy->blocks.begin();
        UELEM(uni10_elem, _package, _type)<uni10_type>*  elem_in;
        UELEM(uni10_elem, _package, _type)<uni10_type>*  elem_out;
        uni10_uint64 Rnum, Cnum;
        Rnum = it_in->second.row_enforce();
        Cnum = it_in->second.col_enforce();
        elem_in = &(it_in->second.elem_enforce());
        elem_out = &(it_out->second.elem_enforce());
        setTranspose(elem_in, &Rnum, &Cnum, elem_out);

      }

    template <typename uni10_type>
      void dagger_nsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para){

        typename std::map<Qnum, Block<uni10_type> >::iterator it_in = src_para->nsy->blocks.begin();
        typename std::map<Qnum, Block<uni10_type> >::iterator it_out = para->nsy->blocks.begin();
        UELEM(uni10_elem, _package, _type)<uni10_type>*  elem_in;
        UELEM(uni10_elem, _package, _type)<uni10_type>*  elem_out;
        uni10_uint64 Rnum, Cnum;
        Rnum = it_in->second.row_enforce();
        Cnum = it_in->second.col_enforce();
        elem_in = &(it_in->second.elem_enforce());
        elem_out = &(it_out->second.elem_enforce());
        setDagger(elem_in, &Rnum, &Cnum, elem_out);

      }

    template <typename uni10_type>
      void conj_nsy(U_para<uni10_type>* para, const U_para<uni10_type>* src_para){

        setConjugate(&src_para->nsy->U_elem, &src_para->nsy->U_elemNum, &para->nsy->U_elem);

      }

    template <typename uni10_type>
      void permute_nsy(const U_para<uni10_type>* T1_para, const std::vector<uni10_int32>& rsp_outin,
          U_para<uni10_type>* T2_para, uni10_bool inorder){

        uni10_int32 bondNum = T1_para->nsy->bonds.size();
        std::vector<uni10_int32> ori_bdDims(bondNum);

        for(uni10_int32 b = 0; b < (uni10_int32)bondNum; b++)
          ori_bdDims[b] = T1_para->nsy->bonds[b].const_getQdegs()[0];

        T1_para->nsy->U_elem.reshape(ori_bdDims, rsp_outin, T2_para->nsy->U_elem, inorder);
        
      }
  
    template <typename uni10_type>
      void addGate_nsy(U_para<uni10_type>* T1_para, const std::vector<_Swap>& swaps){

        return ;

      }

    template <typename uni10_type>
      void traceByRow_nsy(U_para<uni10_type>* Tout_para, const U_para<uni10_type>* Tin_para, uni10_int32 ia, uni10_int32 ib ){

        uni10_int32 bondNum = Tin_para->nsy->bonds.size();
        std::vector<uni10_int32> Q_acc(bondNum, 1);
        for(uni10_int32 b = bondNum - 1; b > 0; b--)
          Q_acc[b - 1] = Q_acc[b] * Tin_para->nsy->bonds[b].const_getQnums().size();

        uni10_int32 tQdim = Tin_para->nsy->bonds[ia].const_getQnums().size();
        uni10_uint64 sB_rDim  = Tin_para->nsy->bonds[ia].dim();
        uni10_uint64 sB_rDim2 = sB_rDim * sB_rDim;
        uni10_error_msg(  !(tQdim == Tin_para->nsy->bonds[ib].const_getQnums().size()), "%s",  "The bonds of the given two labels does not match for trace.");

        Qnum q0(0, PRT_EVEN);
        for(uni10_int32 q = 0; q < tQdim; q++){
          uni10_error_msg(!((Tin_para->nsy->bonds[ia].const_getQnums()[q] * Tin_para->nsy->bonds[ib].const_getQnums()[q] == q0) 
                && (Tin_para->nsy->bonds[ia].const_getQdegs()[q] == Tin_para->nsy->bonds[ib].const_getQdegs()[q]))
              , "%s", "The bonds of the given two labels does not match for trace.");
        }

        typename std::map<Qnum, Block<uni10_type> >::iterator itIn  = Tin_para->nsy->blocks.begin();
        typename std::map<Qnum, Block<uni10_type> >::iterator itOut = Tout_para->nsy->blocks.begin();
       
        uni10_uint64 realElemNum = itOut->second.row() * itOut->second.col();
        for(; itIn != Tin_para->nsy->blocks.end(); itIn++){
          for(uni10_uint64 i = 0; i < realElemNum; i++){
            itOut->second.elem_enforce().__elem[i] = 0;
            for(uni10_uint64 sB_r = 0; sB_r < sB_rDim; sB_r++)
              itOut->second.elem_enforce().__elem[i] += itIn->second.elem_enforce().__elem[i*sB_rDim2+sB_r*sB_rDim+sB_r];
          }
        }

      }

    template <typename uni10_type>
      uni10_type tensorAt_nsy(U_para<uni10_type>* T_para, const uni10_uint64* idxs){

        uni10_int32 bondNum = T_para->nsy->bonds.size();
        std::vector<uni10_int32> Qidxs(bondNum, 0);

        uni10_uint64 bondDim = 1;
        uni10_uint64 idx = 0;
        for(uni10_int32 b = bondNum-1; b >= 0; b--){
          uni10_error_msg(!(idxs[b] < (uni10_uint64)T_para->nsy->bonds[b].dim()), "%s", 
              "The input indices are out of range.");
          idx += bondDim * idxs[b];
          bondDim *= T_para->nsy->bonds[b].dim();

        }

        return T_para->nsy->U_elem[idx];

      }

  };

};

#endif
