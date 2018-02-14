/****************************************************************************
 *  @file Matrix.h
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
 *  @brief Header file for Matrix class
 *  @author Yun-Da Hsieh
 *  @date 2014-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/

#include "uni10/uni10_error.h"
#include "uni10/uni10_api/linalg.h"
#include "uni10/uni10_api/tensor_tools/tensor_tools.h"
#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_permute.h"

namespace uni10{

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(): style(no_sym){

      this->init_para();
      this->meta_link();
      *status = 0;
      this->init();
      (*this->U_elem)[0] = 0;

    };

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(uni10_type val){ 

      this->style = no_sym;
      this->init_para();
      this->meta_link();
      *status = 0;
      this->init();
      (*this->U_elem)[0] = val;
    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const std::vector<Bond>& _bonds, const std::string& _name){

      this->style = check_bonds(_bonds);
      this->init_para();
      this->meta_link();
      *name  = _name;
      *bonds = _bonds;
      *status= 0;
      this->init();

    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const std::vector<Bond>& _bonds, int* _labels, const std::string& _name){

      this->style = check_bonds(_bonds);
      this->init_para();
      this->meta_link();
      *name  = _name;
      *bonds = _bonds;
      *status= 0;
      this->init();
      this->setLabel(_labels);

    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& _labels, const std::string& _name){

      this->style = check_bonds(_bonds);
      this->init_para();
      this->meta_link();
      *name  = _name;
      *bonds = _bonds;
      *status= 0;
      this->init();
      this->setLabel(_labels);

    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const UniTensor& UniT): style(UniT.style){

      if(UniT.paras !=NULL){
        this->init_para();
        this->meta_link();
        this->copy_para(UniT.paras);
        this->initBlocks();
      }else{
        this->init_paras_null();
      }

      ELEMNUM += *this->U_elemNum;
      COUNTER++;
      if(ELEMNUM > MAXELEMNUM)
        MAXELEMNUM = ELEMNUM;
      if(*this->U_elemNum > MAXELEMTEN)
        MAXELEMTEN = *this->U_elemNum;

    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const std::string& fname){

      this->load(fname);

    }

  template <typename uni10_type>
    UniTensor<uni10_type>::UniTensor(const Block<uni10_type>& blk){

      Bond bdi(BD_IN, blk.Rnum);
      Bond bdo(BD_OUT, blk.Cnum);
      this->style = no_sym;
      this->init_para();
      this->meta_link();
      bonds->push_back(bdi);
      bonds->push_back(bdo);
      this->init();
      this->putBlock(blk);

    }

  template <typename uni10_type>
    UniTensor<uni10_type>::~UniTensor(){
      ELEMNUM -= *this->U_elemNum;
      COUNTER--;
      this->free_para();
    }


  template<typename uni10_type>
    UniTensor<uni10_type>& UniTensor<uni10_type>::operator=(UniTensor<uni10_type> const& UniT){

      if(this->paras != NULL)
        this->free_para();

      if(UniT.paras != NULL){
        this->style = UniT.style;
        this->init_para();
        this->meta_link();
        this->copy_para(UniT.paras);
        this->initBlocks();
      }else
        this->init_paras_null();

      return *this; 
    }

  template<typename uni10_type>
    UniTensor<uni10_type>& UniTensor<uni10_type>::operator+=( const UniTensor<uni10_type>& Tb ){
      vectorAdd(this->U_elem, Tb.U_elem, &Tb.U_elem->__elemNum );
      return *this;
    };

  template<typename uni10_type>
    UniTensor<uni10_type>& UniTensor<uni10_type>::operator-=( const UniTensor<uni10_type>& Tb ){
      vectorSub(this->U_elem, Tb.U_elem, &Tb.U_elem->__elemNum );
      return *this;
    };

  template<typename uni10_type>
    UniTensor<uni10_type>& UniTensor<uni10_type>::operator*= (uni10_double64 a){
      vectorScal(&a, this->U_elem, &this->U_elem->__elemNum);
      return *this;
    };

  template <typename uni10_type>
    void UniTensor<uni10_type>::save(const std::string& fname) const{

      FILE* fp = fopen(fname.c_str(), "w");
      uni10_error_msg(!(fp != NULL), "Error in writing to file '%s'", fname.c_str());

      fwrite(&style, 1, sizeof(style), fp);
      uni10_uint64 namelen = name->size();
      fwrite(&namelen, 1, sizeof(namelen), fp);
      fwrite(&(*name)[0], namelen, sizeof((*name)[0]), fp);
      fwrite(&(*status), 1, sizeof(*status), fp);  //OUT: status(4 bytes)
      uni10_uint64 bondNum = bonds->size();
      fwrite(&bondNum, 1, sizeof(bondNum), fp);  //OUT: bondNum(4 bytes)
      uni10_uint64 qnum_sz = sizeof(Qnum);
      fwrite(&qnum_sz, 1, sizeof(qnum_sz), fp);  //OUT: sizeof(Qnum)
      for(uni10_uint64 b = 0; b < bondNum; b++){
        uni10_uint64 num_q = (*bonds)[b].Qnums.size();
        fwrite(&((*bonds)[b].m_type), 1, sizeof(bondType), fp); //OUT: Number of Qnums in the bond(4 bytes)
        fwrite(&num_q, 1, sizeof(num_q), fp);   //OUT: Number of Qnums in the bond(4 bytes)
        fwrite(&((*bonds)[b].Qnums[0]), num_q, qnum_sz, fp);
        fwrite(&((*bonds)[b].Qdegs[0]), num_q, sizeof((*bonds)[b].Qdegs[0]), fp);
      }
      uni10_uint64 num_l = labels->size();
      fwrite(&num_l, 1, sizeof(num_l), fp); //OUT: Number of Labels in the Tensor(4 bytes)
      fwrite(&((*labels)[0]), num_l, sizeof((*labels)[0]), fp);

      if(*status & HAVEELEM)
        U_elem->save(fp);

      fclose(fp);

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::load(const std::string& fname){

      *status = 0;
      this->free_para();
      ELEMNUM -= *this->U_elemNum;
      COUNTER--;

      FILE* fp = fopen(fname.c_str(), "r");
      uni10_error_msg(!(fp != NULL), "Error in opening file '%s'.", fname.c_str());
      uni10_error_msg(!fread(&style, 1, sizeof(style), fp), "%s", "Loading the style of UniTensor is failure. (UniTensor<T>)");
      this->init_para();
      this->meta_link();
      uni10_uint64 namelen;
      uni10_error_msg(!fread(&namelen, 1, sizeof(namelen), fp), "%s", "Loading the namelen of UniTensor is failure. (UniTensor<T>)");
      name->assign(" ", namelen);
      if(namelen != 0)
        uni10_error_msg(!fread(&(*name)[0], namelen, sizeof((*name)[0]), fp), "%s", "Loading the name of UniTensor is failure. (UniTensor<T>)");
      uni10_int st;
      uni10_error_msg(!fread(&st, 1, sizeof(st), fp), "%s", "Loading the status of UniTensor is failure. (UniTensor<T>)");  //OUT: status(4 bytes)
      uni10_uint64 bondNum;
      uni10_error_msg(!fread(&bondNum, 1, sizeof(bondNum), fp), "%s", "Loading the bondNum of UniTensor is failure. (UniTensor<T>)");  //OUT: bondNum(4 bytes)
      uni10_uint64 qnum_sz;
      uni10_error_msg(!fread(&qnum_sz, 1, sizeof(qnum_sz), fp), "%s", "Loading the qnum_sz of UniTensor is failure. (UniTensor<T>)");  //OUT: sizeof(Qnum)
      uni10_error_msg(!(qnum_sz == sizeof(Qnum)), "Error in reading file '%s' in.", fname.c_str());

      for(uni10_uint64 b = 0; b < bondNum; b++){
        uni10_uint64 num_q;
        bondType tp;
        uni10_error_msg(!fread(&tp, 1, sizeof(bondType), fp), "%s", "Loading the bondType of UniTensor is failure. (UniTensor<T>)"); //OUT: Number of Qnums in the bond(4 bytes)
        uni10_error_msg(!fread(&num_q, 1, sizeof(num_q), fp), "%s", "Loading the num_q of UniTensor is failure. (UniTensor<T>)");   //OUT: Number of Qnums in the bond(4 bytes)
        Qnum q0;
        std::vector<Qnum> qnums(num_q, q0);
        uni10_error_msg(!fread(&qnums[0], num_q, qnum_sz, fp), "%s", "Loading the qnums of UniTensor is failure. (UniTensor<T>)");
        std::vector<uni10_int> qdegs(num_q, 0);
        uni10_error_msg(!fread(&(qdegs[0]), num_q, sizeof(qdegs[0]), fp), "%s", "Loading the qdegs of UniTensor is failure. (UniTensor<T>)");
        std::vector<Qnum> tot_qnums;
        for(uni10_uint64 q = 0; q < num_q; q++)
          for(uni10_int d = 0; d < qdegs[q]; d++)
            tot_qnums.push_back(qnums[q]);
        
        Bond bd(tp, tot_qnums);
        bonds->push_back(bd);
      }

      this->init();

      uni10_uint64 num_l;
      uni10_error_msg(!fread(&num_l, 1, sizeof(num_l), fp), "%s", "Loading the num_l of UniTensor is failure. (UniTensor<T>)");  //OUT: Number of Labels in the Tensor(4 bytes)
      uni10_error_msg(!(num_l == bonds->size()), "Error in reading file '%s' in.", fname.c_str());
      labels->assign(num_l, 0);
      if(num_l != 0)
        uni10_error_msg(!fread(&((*labels)[0]), num_l, sizeof((*labels)[0]), fp), "%s", "Loading the labels of UniTensor is failure. (UniTensor<T>)");

      if(st & HAVEELEM){
        U_elem->load(fp);
        *status |= HAVEELEM;
      }

      ELEMNUM += *this->U_elemNum;
      COUNTER++;
      if(ELEMNUM > MAXELEMNUM)
        MAXELEMNUM = ELEMNUM;
      if(*this->U_elemNum > MAXELEMTEN)
        MAXELEMTEN = *this->U_elemNum;

      fclose(fp);

    }

  template <typename uni10_type>
    UniTensor<uni10_type>& UniTensor<uni10_type>::assign(const std::vector<Bond>& _bond ){

      UniTensor<uni10_type> T(_bond);
      *this = T;
      return *this;

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::putBlock(const Block<uni10_type>& mat){
      Qnum q0(0);
      putBlock(q0, mat);
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::putBlock(const Qnum& qnum, const Block<uni10_type>& mat){

      tensor_tools::putBlock(paras, qnum, mat, style);
      *status |= HAVEELEM;

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setName(const std::string& _name){
      *name = _name;
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setLabel(const uni10_int newLabel, const uni10_uint64 idx){
      uni10_error_msg(labels->size() <= idx, "%s", "The bond index is out of the range of vector(labels).");
      (*labels)[idx] = newLabel;
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setLabel(const std::vector<uni10_int>& newLabels){
      uni10_error_msg(!(bonds->size() == newLabels.size()), "%s", "The size of input vector(labels) does not match for the number of bonds.");
      *labels = newLabels;
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setLabel(uni10_int* newLabels){
      std::vector<uni10_int> labels(newLabels, newLabels + bonds->size());
      setLabel(labels);
    }

  template <typename uni10_type>
    const Block<uni10_type>& UniTensor<uni10_type>::const_getBlock()const{
      Qnum q0(0);
      return const_getBlock(q0);
    }

  template <typename uni10_type>
    const Block<uni10_type>& UniTensor<uni10_type>::const_getBlock(const Qnum& qnum)const{
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks->find(qnum);
      if(it == blocks->end()){
        uni10_error_msg(true, "%s", "There is no block with the given quantum number ");
        std::cout << qnum;
      }
      return it->second;
    }

  template <typename uni10_type>
    std::vector<Qnum> UniTensor<uni10_type>::blockQnum()const{
      std::vector<Qnum> keys;
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks->begin();
      for(; it != blocks->end(); it++)
        keys.push_back(it->first);
      return keys;
    }

  template <typename uni10_type>
    Qnum UniTensor<uni10_type>::blockQnum(uni10_uint64 idx)const{

      uni10_error_msg(!(idx < blocks->size()), "Index exceeds the number of the blocks( %ld ).", blocks->size());
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks->begin();
      for(; it != blocks->end(); it++){
        if(idx == 0)
          return it->first;
        idx--;
      }
      return Qnum();
    }

  template <typename uni10_type>
    std::map< Qnum, Matrix<uni10_type> > UniTensor<uni10_type>::getBlocks()const{
      std::map<Qnum, Matrix<uni10_type> > mats;
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks->begin();
      for(; it != blocks->end(); it++){
        Matrix<uni10_type> mat(it->second.Rnum, it->second.Cnum, it->second.elem.__elem);
        mats.insert(std::pair<Qnum, Matrix<uni10_type> >(it->first, mat));
      }
      return std::map< Qnum, Matrix<uni10_type> >();
    }

  template <typename uni10_type>
    Matrix<uni10_type> UniTensor<uni10_type>::getBlock(uni10_bool diag)const{
      Qnum q0(0);
      return getBlock(q0, diag);
    }

  template <typename uni10_type>
    Matrix<uni10_type> UniTensor<uni10_type>::getBlock(const Qnum& qnum, uni10_bool diag)const{
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = blocks->find(qnum);
      if(it == blocks->end()){
        uni10_error_msg(true, "%s", "There is no block with the given quantum number ");
        std::cout<<qnum;
      }
      if(diag)
        return getDiag(it->second);
      else{
        Matrix<uni10_type> mat(it->second.Rnum, it->second.Cnum, it->second.elem.__elem, false);
        return mat;
      }
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setRawElem(const uni10_type* rawElem){

      tensor_tools::setRawElem(paras, rawElem, style);
      *status |= HAVEELEM;

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setRawElem(const std::vector<uni10_type>& rawElem){

      setRawElem(&rawElem[0]);

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setRawElem(const Block<uni10_type>& blk){

      setRawElem(blk.getElem());

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setElem(const uni10_type* _src){
      //UELEM(uni10_elem, _package, _type)<uni10_type> _src(_elem, 1, *U_elemNum, false);
      U_elem->setElem(_src, false);
      *status |= HAVEELEM;
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::setElem(const std::vector<uni10_type>& _elem){
      uni10_error_msg(_elem.size() != *U_elemNum, "%s", "The number of input elements is defferent from the size of the UniTensor");
      setElem(&_elem[0]);
    }

  template <typename uni10_type>
    UniTensor<uni10_type>& UniTensor<uni10_type>::combineBond(const std::vector<uni10_int>& cmbLabels){

      uni10_error_msg((*status & HAVEBOND) == 0, "%s", "There is no bond in the tensor to be combined.");

      if(!(cmbLabels.size() > 1)){
        return *this;
      }

      std::vector<uni10_int> rsp_labels(labels->size(), 0);
      std::vector<uni10_int> reduced_labels(this->labels->size() - cmbLabels.size() + 1, 0);

      std::vector<uni10_int> marked(this->labels->size(), 0);
      std::vector<uni10_int> picked(cmbLabels.size(), 0);
      for(uni10_uint64 p = 0; p < cmbLabels.size(); p++){
        for(uni10_uint64 l = 0; l < this->labels->size(); l++){
          if(cmbLabels[p] == (*this->labels)[l]){
            picked[p] = l;
            marked[l] = 1;
            break;
          }
        }
      }

      uni10_int mark = 0;
      for(uni10_uint64 m = 0; m < marked.size(); m++)
        if(marked[m])
          mark++;

      uni10_error_msg(!(mark == cmbLabels.size()), "%s", "The input labels do not match for the labels of the tensor.");
      uni10_int enc = 0;
      uni10_int enc_r = 0;

      std::vector<Bond> newBonds;
      uni10_int RBnum = 0;
      for(uni10_uint64 l = 0; l < this->labels->size(); l++){
        if(marked[l] && l == (uni10_uint64)picked[0]){
          for(uni10_uint64 ll = 0; ll < cmbLabels.size(); ll++){
            rsp_labels[enc] = cmbLabels[ll];

            enc++;
          }
          std::vector<Bond> tmpBonds;
          for(uni10_uint64 p = 0; p < picked.size(); p++)
            tmpBonds.push_back((*this->bonds)[picked[p]]);
          if((*this->bonds)[picked[0]].type() == BD_IN)
            RBnum += picked.size();
          newBonds.push_back(combine(tmpBonds));
          reduced_labels[enc_r] = (*this->labels)[l];
          enc_r++;
        }
        else if(marked[l] == 0){
          rsp_labels[enc] = (*this->labels)[l];
          reduced_labels[enc_r] = (*this->labels)[l];
          if((*this->bonds)[l].type() == BD_IN)
            RBnum++;
          newBonds.push_back((*this->bonds)[l]);
          enc_r++;
          enc++;
        }
      }

      UniTensor<uni10_type> T_ori;
      permute(T_ori, *this, rsp_labels, RBnum, INPLACE);

      uni10_int isInit = (*status & HAVEELEM);
      this->assign(newBonds);
      this->setLabel(reduced_labels);
      if(isInit){
        this->U_elem->copy(0, *T_ori.U_elem, T_ori.U_elem->__elemNum);
        *status |= HAVEELEM;
      }

      return *this;
    }

  template <typename uni10_type>
    UniTensor<uni10_type>& UniTensor<uni10_type>::combineBond(uni10_int* combined_labels, uni10_int bondNum){

      std::vector<uni10_int> _combined_labels(combined_labels, combined_labels+bondNum);
      return this->combineBond(_combined_labels);

    }

  template<typename uni10_type>
    uni10_type UniTensor<uni10_type>::at(const std::vector<uni10_uint64>& idxs)const{

      uni10_error_msg((*status & HAVEBOND) == 0, "%s", "The tensor is a scalar. Use UniTensor::operator() instead.");
      uni10_error_msg(!(idxs.size() == bonds->size()), "%s", "The size of input indices array does not match with the number of the bonds.");

      return tensor_tools::tensorAt(this->paras, &idxs[0], style);

    }


  template <typename uni10_type>
    void UniTensor<uni10_type>::addGate(const std::vector<_Swap>& swaps){

      tensor_tools::addGate(this->paras, swaps, this->style);

    }

  template <typename uni10_type>
    std::vector<_Swap> UniTensor<uni10_type>::exSwap(const UniTensor<uni10_double64>& Tb)const{

      std::vector<_Swap> swaps;

      if(*status & *Tb.status & HAVEBOND){
        uni10_int bondNumA = labels->size();
        uni10_int bondNumB = Tb.labels->size();
        std::vector<uni10_int> intersect;
        std::vector<uni10_int> left;
        for(uni10_int a = 0; a < bondNumA; a++){
          uni10_bool found = false;
          for(uni10_int b = 0; b < bondNumB; b++)
            if((*labels)[a] == (*Tb.labels)[b])
              found = true;
          if(found)
            intersect.push_back(a);
          else
            left.push_back(a);
        }
        _Swap sp;
        for(uni10_uint64 i = 0; i < intersect.size(); i++)
          for(uni10_uint64 j = 0; j < left.size(); j++){
            sp.b1 = intersect[i];
            sp.b2 = left[j];
            swaps.push_back(sp);
          }
      }
      return swaps;
    }

  template <typename uni10_type>
    std::vector<_Swap> UniTensor<uni10_type>::exSwap(const UniTensor<uni10_complex128>& Tb)const{

      std::vector<_Swap> swaps;

      if(*status & *Tb.status & HAVEBOND){
        uni10_int bondNumA = labels->size();
        uni10_int bondNumB = Tb.labels->size();
        std::vector<uni10_int> intersect;
        std::vector<uni10_int> left;
        for(uni10_int a = 0; a < bondNumA; a++){
          uni10_bool found = false;
          for(uni10_int b = 0; b < bondNumB; b++)
            if((*labels)[a] == (*Tb.labels)[b])
              found = true;
          if(found)
            intersect.push_back(a);
          else
            left.push_back(a);
        }
        _Swap sp;
        for(uni10_uint64 i = 0; i < intersect.size(); i++)
          for(uni10_uint64 j = 0; j < left.size(); j++){
            sp.b1 = intersect[i];
            sp.b2 = left[j];
            swaps.push_back(sp);
          }
      }
      return swaps;
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::printDiagram()const{

      if(this->paras == NULL){
        std::cout<<"This UniTensor has not been initialized."<< std::endl;
        return ;
      }

      if(!(*status & HAVEBOND)){
        if(U_elem->__ongpu)
          std::cout<<"\nScalar: " << U_elem->__elem[0]<<", onGPU";
        else
          std::cout<<"\nScalar: " << U_elem->__elem[0];
        std::cout<<"\n\n";
      }
      else{

        uni10_uint64 row = 0;
        uni10_uint64 col = 0;

        std::vector<Bond>_bonds = *bonds;
        for(uni10_uint64 i = 0; i < _bonds.size(); i++)
          if(_bonds[i].type() == BD_IN)
            row++;
          else
            col++;
        uni10_uint64 layer = std::max(row, col);
        uni10_uint64 nmlen = name->length() + 2;
        uni10_int star = 12 + (14 - nmlen) / 2;
        for(uni10_int s = 0; s < star; s++)
          std::cout << "*";
        if(name->length() > 0)
          std::cout << " " << *name << " ";
        for(uni10_int s = 0; s < star; s++)
          std::cout <<"*";
        std::cout<<std::endl;

        if(U_elem->__uni10_typeid == 1)
          std::cout << "REAL" << std::endl;
        else if(U_elem->__uni10_typeid == 2)
          std::cout << "COMPLEX" << std::endl;

        if(U_elem->__ongpu)
          std::cout<<"\n                 onGPU";
        std::cout << "\n             ____________\n";
        std::cout << "            |            |\n";
        uni10_uint64 llab = 0;
        uni10_uint64 rlab = 0;
        char buf[128];
        for(uni10_uint64 l = 0; l < layer; l++){
          if(l < row && l < col){
            llab = (*labels)[l];
            rlab = (*labels)[row + l];
            sprintf(buf, "    %5ld___|%-4d    %4d|___%-5ld\n", llab, _bonds[l].dim(), _bonds[row + l].dim(), rlab);
            std::cout<<buf;
          }
          else if(l < row){
            llab = (*labels)[l];
            sprintf(buf, "    %5ld___|%-4d    %4s|\n", llab, _bonds[l].dim(), "");
            std::cout<<buf;
          }
          else if(l < col){
            rlab = (*labels)[row + l];
            sprintf(buf, "    %5s   |%4s    %4d|___%-5ld\n", "", "", _bonds[row + l].dim(), rlab);
            std::cout << buf;
          }
          std::cout << "            |            |   \n";
        }
        std::cout << "            |____________|\n";

        std::cout << "\n================BONDS===============\n";
        for(uni10_uint64 b = 0; b < _bonds.size(); b++)
          std::cout << _bonds[b];

        std::cout << "\n\nTotal elemNum: "<<(*U_elemNum)<<std::endl;
        std::cout << "====================================\n";

      }
    }

  template<typename uni10_type>
    Matrix<uni10_type> UniTensor<uni10_type>::getRawElem()const{

      if(*status & HAVEBOND && *status & HAVEELEM){
        uni10_int bondNum = bonds->size();
        uni10_uint64 rowNum = 1;
        uni10_uint64 colNum = 1;
        for(std::vector<Bond>::const_iterator it = bonds->begin(); it != bonds->end(); ++it){
          if(it->type() == BD_IN)
            rowNum *= it->dim();
          else
            colNum *= it->dim();
        }
        std::vector<uni10_uint64> idxs(bondNum, 0);
        uni10_int bend;
        std::vector<uni10_type> rawElem;
        while(1){
          rawElem.push_back(at(idxs));
          //std::cout << "testing: " << at(idxs) << std::endl;
          for(bend = bondNum - 1; bend >= 0; bend--){
            idxs[bend]++;
            if(idxs[bend] < (uni10_uint64)(*bonds)[bend].dim())
              break;
            else
              idxs[bend] = 0;
          }
          if(bend < 0)
            break;
        }
        return Matrix<uni10_type>(rowNum, colNum, &rawElem[0]);
      }
      else if(*status & HAVEELEM){
        return Matrix<uni10_type>(1, 1, this->U_elem->__elem);
      }

      return Matrix<uni10_type>();

    }

  template<typename uni10_type>
    std::string UniTensor<uni10_type>::printRawElem(bool print)const{

      std::ostringstream os;
      if(*status & HAVEBOND && *status & HAVEELEM){
        uni10_int bondNum = bonds->size();
        std::vector<Bond> ins;
        std::vector<Bond> outs;
        for(std::vector<Bond>::const_iterator it = bonds->begin(); it != bonds->end(); ++it){
          if(it->type() == BD_IN)
            ins.push_back(*it);
          else
            outs.push_back(*it);
        }
        if(ins.size() == 0 || outs.size() == 0)
          os<<getRawElem();
        else{
          Bond rBond = combine(ins);
          Bond cBond = combine(outs);
          std::vector<Qnum> rowQ = rBond.Qlist();
          std::vector<Qnum> colQ = cBond.Qlist();
          uni10_uint64 colNum = cBond.dim();
          std::vector<uni10_uint64> idxs(bondNum, 0);

          os<< "     ";
          for(uni10_uint64 q = 0; q < colQ.size(); q++)
            os<< "   " << std::setw(2) << colQ[q].U1() << "," << colQ[q].prt();
          os<< std::endl << std::setw(5) << "" << std::setw(colQ.size() * 7 + 2) <<std::setfill('-')<<"";
          os<<std::setfill(' ');
          uni10_int cnt = 0;
          uni10_int r = 0;
          uni10_int bend;
          while(1){
            if(cnt % colNum == 0){
              os<<"\n    |\n" << std::setw(2) << rowQ[r].U1() << "," << rowQ[r].prt() << "|";
              r++;
            }
            os<< std::setw(7) << std::fixed << std::setprecision(3) << at(idxs);
            for(bend = bondNum - 1; bend >= 0; bend--){
              idxs[bend]++;
              if(idxs[bend] < (uni10_uint64)(*bonds)[bend].dim())
                break;
              else
                idxs[bend] = 0;
            }
            cnt++;
            if(bend < 0)
              break;
          }
          os <<"\n    |\n";
        }
      }
      if(*status & HAVEELEM){
        //os<<"\nScalar: " << c_elem[0]<<"\n\n";
      }
      else
        os<<"NO ELEMENT IN THE TENSOR!!!\n";

      if(print){
        std::cout<<os.str();
        return "";
      }

      return os.str();

    }

  template<typename uni10_type>
    std::string UniTensor<uni10_type>::profile(bool print){

      std::ostringstream os;
      os<<"\n===== Tensor profile =====\n";
      os<<"Existing Tensors: " << COUNTER << std::endl;
      os<<"Allocated Elements: " << ELEMNUM << std::endl;
      os<<"Max Allocated Elements: " << MAXELEMNUM << std::endl;
      os<<"Max Allocated Elements for a Tensor: " << MAXELEMTEN << std::endl;
      os<<"============================\n\n";
      if(print){
        std::cout<<os.str();
        return "";
      }
      return os.str();

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::init_para(){

      paras = tensor_tools::init_para(paras, style);

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::copy_para(U_para<uni10_type>* src_para){

      tensor_tools::copy_para(paras, src_para, style);

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::meta_link(){

      if(style == 0){

        this->name      = &paras->nsy->name;
        this->bonds     = &paras->nsy->bonds;
        this->labels    = &paras->nsy->labels;
        this->RBondNum  = &paras->nsy->RBondNum;  
        this->RQdim     = &paras->nsy->RQdim;
        this->CQdim     = &paras->nsy->CQdim;
        this->U_elemNum = &paras->nsy->U_elemNum;    
        this->blocks    = &paras->nsy->blocks;
        this->U_elem    = &paras->nsy->U_elem;
        this->status    = &paras->nsy->status;

      }
      else if(style == 1){

        this->name      = &paras->bsy->name;
        this->bonds     = &paras->bsy->bonds;
        this->labels    = &paras->bsy->labels;
        this->RBondNum  = &paras->bsy->RBondNum;  
        this->RQdim     = &paras->bsy->RQdim;
        this->CQdim     = &paras->bsy->CQdim;
        this->U_elemNum = &paras->bsy->U_elemNum;    
        this->blocks    = &paras->bsy->blocks;
        this->U_elem    = &paras->bsy->U_elem;
        this->status    = &paras->bsy->status;

      }
      else if(style == 2){
        //name = &paras->ssy->name;
        uni10_error_msg(true, "%s", "Developping!!!");
      }

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::init(){
      // You should init_para() first. Then you can use this function to initialize the UniTensor.
      tensor_tools::init(paras, style);
      this->check_status();

      ELEMNUM += *this->U_elemNum;
      COUNTER++;
      if(ELEMNUM > MAXELEMNUM)
        MAXELEMNUM = ELEMNUM;
      if(*this->U_elemNum > MAXELEMTEN)
        MAXELEMTEN = *this->U_elemNum;

    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::check_status(){

      if(bonds->size())
          *status |= HAVEBOND;
      else
          *status |= HAVEELEM;
      
    }

  template <typename uni10_type>
    void UniTensor<uni10_type>::initBlocks(){

      tensor_tools::initBlocks(paras, style);

    }


  template <typename uni10_type>
    void UniTensor<uni10_type>::free_para(){
      // You should init_para() first. Then you can use this function to initialize the UniTensor.
      tensor_tools::free_para(paras, style);
    }

  template<typename uni10_type> 
    void UniTensor<uni10_type>::init_paras_null(){

      this->paras = NULL;
      this->name = NULL;
      this->bonds = NULL;
      this->labels = NULL;
      this->RBondNum = NULL;     
      this->RQdim = NULL;
      this->CQdim = NULL;
      this->U_elemNum = NULL;
      this->blocks = NULL;
      this->U_elem = NULL;
      this->status = NULL;    


    };

  template <typename uni10_type> 
    contain_type UniTensor<uni10_type>::check_bonds(const std::vector<Bond>& _bonds)const{

      contain_type s;
      uni10_bool withoutSymmetry = true;
      for(uni10_uint64 b = 0; b < _bonds.size(); b++){
        if(_bonds[b].const_getQnums().size() != 1){
          withoutSymmetry = false;
          break;
        }
      }
      // spar_sym are developping
      s = withoutSymmetry ? no_sym : blk_sym;
      return s;

    }

  template<typename uni10_type> 
    void UniTensor<uni10_type>::set_zeros(){

      this->U_elem->set_zeros();
      *status |= HAVEELEM;

    };

  template<typename uni10_type> 
    void UniTensor<uni10_type>::set_zeros(const Qnum& qnum){

      typename std::map<Qnum, Block<uni10_type> >::iterator it = this->blocks->find(qnum);

      uni10_error_msg(it == this->blocks->end(), "%s", "There is no block with the given quantum number.");

      Block<uni10_type>& block = it->second;
      block.elem.set_zeros();
      *status |= HAVEELEM;

    };

  template<typename uni10_type> 
    void UniTensor<uni10_type>::identity(){

      typename std::map<Qnum, Block<uni10_type> >::iterator it = this->blocks->begin();
      for(; it != this->blocks->end(); it++)
        setIdentity(&it->second.elem, &it->second.diag, &it->second.Rnum, &it->second.Cnum);
      *status |= HAVEELEM;

    };

  template<typename uni10_type> 
    void UniTensor<uni10_type>::identity(const Qnum& qnum){

      typename std::map<Qnum, Block<uni10_type> >::iterator it = this->blocks->find(qnum);
      uni10_error_msg(it == this->blocks->end(), "%s", "There is no block with the given quantum number.");
      Block<uni10_type>& block = it->second;
      setIdentity(&block.elem, &block.diag, &block.Rnum, &block.Cnum);
      *status |= HAVEELEM;

    };

  template<typename uni10_type> 
    void UniTensor<uni10_type>::randomize(char UorN, uni10_double64 dn_mu, uni10_double64 up_var, uni10_int64 seed){

      typename std::map<Qnum, Block<uni10_type> >::iterator it = this->blocks->begin();
      if(UorN == 'U'){
        uni10_double64 dn = std::min(dn_mu, up_var);
        uni10_double64 up = std::max(dn_mu, up_var);
        for(; it != this->blocks->end(); it++)
          setUniformRand(&it->second.elem, &it->second.diag, &it->second.Rnum, &it->second.Cnum, &dn, &up, &seed);
      }else if(UorN == 'N'){
        for(; it != this->blocks->end(); it++)
          setNormalRand(&it->second.elem, &it->second.diag, &it->second.Rnum, &it->second.Cnum, &dn_mu, &up_var, &seed);
      }else
        uni10_error_msg(true, "%s", "Wrong flag. Hint: The fisrt parameter must be 'N' or 'U'");
      *status |= HAVEELEM;

    };

  template<typename uni10_type> 
    void UniTensor<uni10_type>::randomize(const Qnum& qnum, char UorN, uni10_double64 dn_mu, uni10_double64 up_var, uni10_int64 seed){

      typename std::map<Qnum, Block<uni10_type> >::iterator it = this->blocks->find(qnum);
      uni10_error_msg(it == this->blocks->end(), "%s", "There is no block with the given quantum number.");
      Block<uni10_type>& block = it->second;
      if(UorN == 'U'){
        uni10_double64 dn = std::min(dn_mu, up_var);
        uni10_double64 up = std::max(dn_mu, up_var);
        setUniformRand(&block.elem, &block.diag, &block.Rnum, &block.Cnum, &dn, &up, &seed);
      }
      else if(UorN == 'N')
        setNormalRand(&block.elem, &block.diag, &block.Rnum, &block.Cnum, &dn_mu, &up_var, &seed);
      else
        uni10_error_msg(true, "%s", "Wrong flag. Hint: The fisrt parameter must be 'N' or 'U'");
      *status |= HAVEELEM;

    };

  template<typename uni10_type> 
    void UniTensor<uni10_type>::orthoRand(char UorN, uni10_double64 dn_mu, uni10_double64 up_var, uni10_int64 seed){

      typename std::map<Qnum, Block<uni10_type> >::iterator it = this->blocks->begin();
      Matrix<uni10_type> U, S, vT;
      Matrix<uni10_type>* null_mat = NULL;
      if(UorN == 'U'){
        uni10_double64 dn = std::min(dn_mu, up_var);
        uni10_double64 up = std::max(dn_mu, up_var);
        for(; it != this->blocks->end(); it++){
          setUniformRand(&it->second.elem, &it->second.diag, &it->second.Rnum, &it->second.Cnum, &dn, &up, &seed);
          if(it->second.Rnum < it->second.Cnum){
            svd(it->second, *null_mat, S, vT, INPLACE);
            this->putBlock(it->first, vT);
          }else{
            svd(it->second, U, S, *null_mat, INPLACE);
            this->putBlock(it->first, U);
          }
        }
      }else if(UorN == 'N'){
        for(; it != this->blocks->end(); it++){
          setNormalRand(&it->second.elem, &it->second.diag, &it->second.Rnum, &it->second.Cnum, &dn_mu, &up_var, &seed);
          if(it->second.Rnum < it->second.Cnum){
            svd(it->second, *null_mat, S, vT, INPLACE);
            this->putBlock(it->first, vT);
          }else{
            svd(it->second, U, S, *null_mat, INPLACE);
            this->putBlock(it->first, U);
          }
        }
      }else
        uni10_error_msg(true, "%s", "Wrong flag. Hint: The fisrt parameter must be 'N' or 'U'");
      *status |= HAVEELEM;

    };

  template<typename uni10_type> 
    void UniTensor<uni10_type>::orthoRand(const Qnum& qnum, char UorN, uni10_double64 dn_mu, uni10_double64 up_var, uni10_int64 seed){

      typename std::map<Qnum, Block<uni10_type> >::iterator it = this->blocks->find(qnum);
      uni10_error_msg(it == this->blocks->end(), "%s", "There is no block with the given quantum number.");
      Matrix<uni10_type> U, S, vT;
      Matrix<uni10_type>* null_mat = NULL;
      Block<uni10_type>& block = it->second;
      if(UorN == 'U'){
        uni10_double64 dn = std::min(dn_mu, up_var);
        uni10_double64 up = std::max(dn_mu, up_var);
        setUniformRand(&block.elem, &block.diag, &block.Rnum, &block.Cnum, &dn, &up, &seed);
        if(it->second.Rnum < it->second.Cnum){
          svd(it->second, *null_mat, S, vT, INPLACE);
          this->putBlock(it->first, vT);
        }else{
          svd(it->second, U, S, *null_mat, INPLACE);
          this->putBlock(it->first, U);
        }
      }
      else if(UorN == 'N'){
        setNormalRand(&block.elem, &block.diag, &block.Rnum, &block.Cnum, &dn_mu, &up_var, &seed);
        if(it->second.Rnum < it->second.Cnum){
          svd(it->second, *null_mat, S, vT, INPLACE);
          this->putBlock(it->first, vT);
        }else{
          svd(it->second, U, S, *null_mat, INPLACE);
          this->putBlock(it->first, U);
        }
      }
      else
        uni10_error_msg(true, "%s", "Wrong flag. Hint: The fisrt parameter must be 'N' or 'U'");
      *status |= HAVEELEM;

    };

  template<typename uni10_type> 
    void UniTensor<uni10_type>::clear(){
      *this = UniTensor();
      (*status) &= ~HAVEELEM;
    }

  UniTensor<uni10_complex128>& operator*=(UniTensor<uni10_complex128>& t1, uni10_complex128 a){
      vectorScal(&a, t1.U_elem, &t1.U_elem->__elemNum);
      return t1;
  }

  template class UniTensor<uni10_double64>;
  template class UniTensor<uni10_complex128>;
}

