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
#include <math.h>

#include "uni10/uni10_error.h"
#include "uni10/uni10_api/Matrix.h"

namespace uni10{


  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(): Block<uni10_type>(){};

  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _diag): Block<uni10_type>(_Rnum, _Cnum, _diag){
      init();
    };

  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_type* _src, uni10_bool _diag): Block<uni10_type>(_Rnum, _Cnum, _diag){
      init(_src);
    };

  // Copy constructor
  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(Matrix const& _m): Block<uni10_type>(_m.Rnum, _m.Cnum, _m.diag){
      this->init(_m.elem);
    };

  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(Block<uni10_type> const& _b): Block<uni10_type>(_b.Rnum, _b.Cnum, _b.diag){
      this->init( _b.elem);
    };

  template <typename uni10_type>
    Matrix<uni10_type>::Matrix(const std::string& fname){

      FILE* fp = fopen(fname.c_str(), "r");
      uni10_error_msg(!fread(&this->Rnum, sizeof(this->Rnum), 1, fp), "%s", "Loading Rnum is failure. (Matrix<T>)");
      uni10_error_msg(!fread(&this->Cnum, sizeof(this->Cnum), 1, fp), "%s", "Loading Cnum is failure. (Matrix<T>)");
      uni10_error_msg(!fread(&this->diag, sizeof(this->diag), 1, fp), "%s", "Loading diag is failure. (Matrix<T>)");
      this->elem.load(fp);
      fclose(fp);

    };

  template <typename uni10_type>
    Matrix<uni10_type>::~Matrix(){};

  template<typename uni10_type>
    Matrix<uni10_type>& Matrix<uni10_type>::operator=(const Matrix<uni10_type>& _m){
      this->Rnum = _m.Rnum;
      this->Cnum = _m.Cnum;
      this->diag = _m.diag;
      this->init(_m.elem);
      return *this;
    }

  template<typename uni10_type>
    Matrix<uni10_type>& Matrix<uni10_type>::operator=(const Block<uni10_type>& _b){
      this->Rnum = _b.Rnum;
      this->Cnum = _b.Cnum;
      this->diag = _b.diag;
      this->init(_b.elem);
      return *this;
    }

  template<typename uni10_type>
    Matrix<uni10_type>& Matrix<uni10_type>::operator+=(const Matrix<uni10_type>& _m){
      matrixAdd(&this->elem, &this->diag, &_m.elem, &_m.diag, &_m.Rnum, &_m.Cnum );
      this->diag = (this->diag && _m.diag);
      return *this;
    }

  template<typename uni10_type>
    Matrix<uni10_type>& Matrix<uni10_type>::operator-=(const Matrix<uni10_type>& _m){
      matrixSub(&this->elem, &this->diag, &_m.elem, &_m.diag, &_m.Rnum, &_m.Cnum );
      this->diag = (this->diag && _m.diag);
      return *this;
    };

  template<typename uni10_type>
    Matrix<uni10_type>& Matrix<uni10_type>::operator*=(uni10_double64 a){ 
      vectorScal(&a, &this->elem, &this->elem.__elemNum);
      return *this;
    }

  template<typename uni10_type>
    Matrix<uni10_type>& Matrix<uni10_type>::operator*=(const Matrix<uni10_type>& _m){            
      matrixMul(&this->elem, &this->diag, &_m.elem, &_m.diag, &_m.Rnum, &_m.Cnum );
      this->diag = (this->diag || _m.diag);
      return *this;
    }

  template <typename uni10_type>
    void Matrix<uni10_type>::assign(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag){
      this->diag = _isdiag;
      this->Rnum = _Rnum;
      this->Cnum = _Cnum;
      this->init();
    };

  template <typename uni10_type>
    void Matrix<uni10_type>::load(const std::string& fname){

      FILE* fp = fopen(fname.c_str(), "r");
      uni10_error_msg(!fread(&this->Rnum, sizeof(this->Rnum), 1, fp), "%s", "Loading Rnum is failure. (Matrix<T>)");
      uni10_error_msg(!fread(&this->Cnum, sizeof(this->Cnum), 1, fp), "%s", "Loading Cnum is failure. (Matrix<T>)");
      uni10_error_msg(!fread(&this->diag, sizeof(this->diag), 1, fp), "%s", "Loading diag is failure. (Matrix<T>)");
      this->elem.load(fp);
      fclose(fp);

    };

  template <typename uni10_type>

    void Matrix<uni10_type>::setElem(const uni10_type* src, bool src_ongpu){

      uni10_error_msg( src_ongpu, "%s", " The source pointer is on the device. Please install the MAGMA or CUDAONLY gpu version instead.");

      if(this->elem.__uni10_typeid != UNI10_TYPE_ID(uni10_type)){

        uni10_error_msg( true, "%s", " Developping !!!");

      }

      this->elem.setElem(src);

    };

  template <typename uni10_type>
    void Matrix<uni10_type>::setElem(const std::vector<uni10_type>& elem, bool src_ongpu){

      uni10_error_msg( src_ongpu, "%s", " The source pointer is on the device. Please install the MAGMA or CUDAONLY gpu version instead.");

      uni10_error_msg(this->diag == false && this->Rnum*this->Cnum != elem.size(),
          "Number of the input elements is: %ld, and it doesn't match to the size of matrix: %ld", elem.size(), this->Rnum*this->Cnum);

      uni10_error_msg(this->diag == true && fmin(this->Rnum, this->Cnum) != elem.size(), 
          "Number of the input elements is: %ld, and it doesn't match to the size of matrix: %.0f", elem.size(), fmin(this->Rnum,this->Cnum));

      setElem(&elem[0], src_ongpu);

    };

  template <typename uni10_type>
    void Matrix<uni10_type>::uni10_elem_free(){
      uni10_error_msg(true, "%s","Developping");
    };

  template <typename uni10_type>
    void Matrix<uni10_type>::set_elem_null(){
      uni10_error_msg(true, "%s", "Developping");
    };

  template <typename uni10_type>
    void Matrix<uni10_type>::init(const uni10_type* elem){

        this->elem.init(this->Rnum, this->Cnum, this->diag, elem);

    };

  template <typename uni10_type>
    void Matrix<uni10_type>::init(const UELEM(uni10_elem, _package, _type)<uni10_type>& _m){

        this->elem.init(this->Rnum, this->Cnum, this->diag,_m);

    };

  template <typename uni10_type>
    void Matrix<uni10_type>::setDiag(const uni10_bool _isdiag){

        this->diag = _isdiag;

    };

  template <typename uni10_type>
    void Matrix<uni10_type>::set_zeros(){

      this->elem.set_zeros();

    };

  template <typename uni10_type>
    void Matrix<uni10_type>::identity(){

      uni10_error_msg(this->elemNum() == 0, "%s", "The matrix has not been initialized!!!");

      this->diag = true;
      if(!this->empty())
        this->elem.clear();
      this->elem.init(this->Rnum, this->Cnum, this->diag);

      for(uni10_uint64 i = 0; i < this->elem.__elemNum; i++)
        this->elem.__elem[i] = 1.;

    };
  
  template<typename uni10_type> 
    void Matrix<uni10_type>::randomize(char UorN, uni10_double64 dn_mu, uni10_double64 up_var, uni10_int64 seed){

      if(UorN == 'U'){
        uni10_double64 dn = std::min(dn_mu, up_var);
        uni10_double64 up = std::max(dn_mu, up_var);
        setUniformRand(&this->elem, &this->diag, &this->Rnum, &this->Cnum, &dn, &up, &seed);
      }else if(UorN == 'N'){
          setNormalRand(&this->elem, &this->diag, &this->Rnum, &this->Cnum, &dn_mu, &up_var, &seed);
      }else
        uni10_error_msg(true, "%s", "Wrong flag. Hint: The fisrt parameter must be 'N' or 'U'");

    };

  template<typename uni10_type> 
    void Matrix<uni10_type>::orthoRand(char UorN, uni10_double64 dn_mu, uni10_double64 up_var, uni10_int64 seed){

      if(UorN == 'U'){
        uni10_double64 dn = std::min(dn_mu, up_var);
        uni10_double64 up = std::max(dn_mu, up_var);
        setUniformRand(&this->elem, &this->diag, &this->Rnum, &this->Cnum, &dn, &up, &seed);
      }
      else if(UorN == 'N'){
        setNormalRand(&this->elem, &this->diag, &this->Rnum, &this->Cnum, &dn_mu, &up_var, &seed);
        Matrix<uni10_type> tmp = *this;
      }
      else
        uni10_error_msg(true, "%s", "Wrong flag. Hint: The fisrt parameter must be 'N' or 'U'");

      uni10_uint64 min = std::min(this->Rnum, this->Cnum);
      UELEM(uni10_elem, _package, _type)<uni10_type>* U_elem  = NULL;     // pointer to a real matrix
      UELEM(uni10_elem, _package, _type)<uni10_type>* vT_elem = NULL;     // pointer to a real matrix
      UELEM(uni10_elem, _package, _type)<uni10_type> S(this->Rnum, this->Cnum, true);
      if(this->Rnum < this->Cnum){
        vT_elem = new UELEM(uni10_elem, _package, _type)<uni10_type>(min, this->Cnum, this->diag);
        matrixSVD(&this->elem, &this->diag, &this->Rnum, &this->Cnum, U_elem, &S, vT_elem);
        this->elem.copy(*vT_elem);
      }else{
        U_elem = new UELEM(uni10_elem, _package, _type)<uni10_type>(this->Rnum, min, this->diag);
        matrixSVD(&this->elem, &this->diag, &this->Rnum, &this->Cnum, U_elem, &S, vT_elem);
        this->elem.copy(*U_elem);
      }

      if(U_elem != NULL)
        delete U_elem;
      else if(vT_elem != NULL)
        delete vT_elem;

    };

  Matrix<uni10_complex128>& operator*=(Matrix<uni10_complex128>& _m, uni10_complex128 a){ 
    vectorScal(&a, &_m.elem, &_m.elem.__elemNum);
    return _m;
  }

  template class Matrix<uni10_double64>;
  template class Matrix<uni10_complex128>;

}  /* namespace uni10 */
