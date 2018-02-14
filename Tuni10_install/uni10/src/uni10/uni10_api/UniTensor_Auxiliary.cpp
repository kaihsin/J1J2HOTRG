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
 *  @author Yun-Hsuan Chou
 *  @date 2014-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/

#include "uni10/uni10_error.h"
#include "uni10/uni10_api/linalg.h"
#include "uni10/uni10_api/tensor_tools/tensor_tools.h"
#include "uni10/uni10_api/uni10_hirnk_linalg_inplace/uni10_hirnk_linalg_inplace_permute.h"

namespace uni10{

  template<> template<>
    UniTensor<uni10_complex128>::UniTensor(UniTensor<uni10_double64> const& UniT): style(UniT.style){

      if(UniT.paras != NULL){

        this->init_para();
        this->meta_link();
        *name  = *UniT.name;
        *bonds = *UniT.bonds;
        *status= 0;
        this->init();
        this->setLabel(*UniT.labels);
        this->U_elem->copy(*(UniT.U_elem));
        *status |= HAVEELEM;

      }else{
        this->init_paras_null();
      }

    }

  template<> template<>
    UniTensor<uni10_double64>::UniTensor(UniTensor<uni10_complex128> const& UniT): style(UniT.style){

      if(UniT.paras != NULL){
        this->init_para();
        this->meta_link();
        *name  = *UniT.name;
        *bonds = *UniT.bonds;
        *status= 0;
        this->init();
        this->setLabel(*UniT.labels);
        this->U_elem->copy(*(UniT.U_elem));
        *status |= HAVEELEM;
      }
      else{
        this->init_paras_null();
      }

    }

  template<> template<>
    UniTensor<uni10_complex128>& UniTensor<uni10_complex128>::operator=(UniTensor<uni10_double64> const& UniT){

      if(this->paras != NULL)
        this->free_para();

      if(UniT.paras != NULL){
        this->style = UniT.style;
        this->init_para();
        this->meta_link();
        *name  = *UniT.name;
        *bonds = *UniT.bonds;
        *status= 0;
        this->init();
        this->setLabel(*UniT.labels);
        this->U_elem->copy(*(UniT.U_elem));
        *status |= HAVEELEM;
      }else{
        this->init_paras_null();
      }

      return *this; 

    }

  template<> template<>
    UniTensor<uni10_double64>& UniTensor<uni10_double64>::operator=(UniTensor<uni10_complex128> const& UniT){

      if(this->paras != NULL)
        this->free_para();

      if(UniT.paras != NULL){
        this->style = UniT.style;
        this->init_para();
        this->meta_link();
        *name  = *UniT.name;
        *bonds = *UniT.bonds;
        *status= 0;
        this->init();
        this->setLabel(*UniT.labels);
        this->U_elem->copy(*(UniT.U_elem));
        *status |= HAVEELEM;
      }else{
        this->init_paras_null();
      }
      return *this; 

    }

  UniTensor<uni10_complex128>& operator+=(UniTensor<uni10_complex128>& t1, const UniTensor<uni10_double64>& t2){

    vectorAdd(t1.U_elem, t2.U_elem, &t1.U_elem->__elemNum);
    return t1;

  }
  
  UniTensor<uni10_complex128>& operator-=(UniTensor<uni10_complex128>& t1, const UniTensor<uni10_double64>& t2){

    vectorSub(t1.U_elem, t2.U_elem, &t1.U_elem->__elemNum);
    return t1;

  }

};
