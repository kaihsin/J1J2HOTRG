/****************************************************************************
 *  @file Block.cpp
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
 *  @brief Implementation file of Block class
 *  @author Yun-Da Hsieh
 *  @date 2014-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/
#include "uni10/uni10_api/Block.h"


namespace uni10{

  template<typename uni10_type>
    Block<uni10_type>::Block(): Rnum(0), Cnum(0), diag(false){};

  template<typename uni10_type>
    Block<uni10_type>::Block(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _diag): Rnum(_Rnum), Cnum(_Cnum), diag(_diag){

    }

  template<typename uni10_type>
    Block<uni10_type>::Block(const Block& _b): Rnum(_b.Rnum), Cnum(_b.Cnum), diag(_b.diag){
      elem.__uni10_typeid = _b.elem.__uni10_typeid;
      elem.__ongpu = _b.elem.__ongpu;
      elem.__elemNum = _b.elem.__elemNum;
      elem.__elem = _b.elem.__elem;
    }

  template<typename uni10_type>
    Block<uni10_type>::~Block(){}

  template<typename uni10_type>
    uni10_uint64 Block<uni10_type>::row()const{return Rnum;}

  template<typename uni10_type>
    uni10_uint64 Block<uni10_type>::col()const{return Cnum;}

  template<typename uni10_type>
    bool Block<uni10_type>::isDiag()const{return diag;}

  template<typename uni10_type>
    uni10_uint64 Block<uni10_type>::elemNum()const{ return elem.__elemNum; }

  template<typename uni10_type>
    int Block<uni10_type>::typeID()const{ return elem.__uni10_typeid;}

  template<typename uni10_type>
    uni10_type* Block<uni10_type>::getElem()const{ return elem.__elem; }

  template<typename uni10_type>
    bool Block<uni10_type>::empty()const{
      
      return this->elem.empty();

    }

  template<typename uni10_type>
    void Block<uni10_type>::save(const std::string& fname)const{

      FILE* fp = fopen(fname.c_str(), "w");
      fwrite(&Rnum, sizeof(Rnum), 1, fp);
      fwrite(&Cnum, sizeof(Cnum), 1, fp);
      fwrite(&diag, sizeof(diag), 1, fp);
      elem.save(fp);
      fclose(fp);

    }

  template<typename uni10_type>
    uni10_type Block<uni10_type>::at(uni10_uint64 r, uni10_uint64 c)const{
      if(diag){
        if(r == c)
          return elem.__elem[c];
        else
          return 0.;
      }
      else
        return elem.__elem[r*Cnum+c];
    }

#if defined(GPU)
  template<typename uni10_type>
    bool Block<uni10_type>::isOngpu()const{return elem.ongpu;}
#endif

  template class Block<uni10_double64>;
  template class Block<uni10_complex128>;

};
