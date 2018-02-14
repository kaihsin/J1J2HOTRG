/****************************************************************************
*  @file Block.h
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
*  @brief Header file for Block class
*  @author Yun-Da Hsieh
*  @author Yun-Hsuan Chou
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef __UNI10_BLOCK_H__
#define __UNI10_BLOCK_H__

#include <assert.h>

#include <iostream>
#include <iomanip>
#include <vector>

#include "uni10/uni10_type.h"
#include "uni10/uni10_elem.h"
#include "uni10/uni10_elem_linalg.h"

namespace uni10{

  enum UNI10_INPLACE{
    INPLACE = 1        //< Modified elements in the Matrix directly.
  };

  template<typename uni10_type>
    class UniTensor;

  template<typename uni10_type>
    class Matrix;

  template<typename uni10_type>
    class Block;

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, const Block<uni10_type>& _b);

  template<typename T, typename U>
    uni10_bool operator==(const Block<T>& m1, const Block<U>& m2);

  template<typename T, typename U>
    uni10_bool operator!=(const Block<T>& m1, const Block<U>& m2);

  /// @class Block
  /// @brief Base class for Matrix.
  ///
  /// A Block holds a reference to a Matrix. The Block constructor does not allocate memory. Memory allocation
  /// should be done through Matrix.
  ///
  /// @see \ref Matrix, UniTensor

  template<typename uni10_type>
    class Block{  

      protected:

        UELEM(uni10_elem, _package, _type)<uni10_type> elem;     // pointer to a real matrix

        uni10_uint64 Rnum;

        uni10_uint64 Cnum;

        uni10_bool diag;

      public:

        // Four except enforce functions which are designed for tensor_tools.
        uni10_uint64& row_enforce(){return Rnum;};

        uni10_uint64& col_enforce(){return Cnum;};

        uni10_bool& diag_enforce(){return diag;};
        
        UELEM(uni10_elem, _package, _type)<uni10_type>& elem_enforce(){return elem;};

        const UELEM(uni10_elem, _package, _type)<uni10_type>& const_elem_enforce()const{return elem;};

        explicit Block();

        explicit Block(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _diag = false);

        explicit Block(const Block& _b);

        virtual ~Block();

        uni10_uint64 row()const;

        uni10_uint64 col()const;

        bool isDiag()const;

        void save(const std::string& fname)const;

        bool empty()const;              // --> uni10_elem().empty()

        uni10_uint64 elemNum()const;    // --> uni10_elem().elemNum()

        uni10_int32 typeID()const;      // --> uni10_elem().typeid()

        bool isOngpu()const;            // --> uni10_elem().isOngpu()

        uni10_type* getElem()const;     // --> uni10_elem().getElem()

        uni10_type at(uni10_uint64 i, uni10_uint64 j)const;

        Block& operator=(const Block& _m){
          //std::cout << "BLock QQ =  !!\n\n";
          this->Rnum = _m.Rnum;
          this->Cnum = _m.Cnum;
          this->diag = _m.diag;
          elem.__uni10_typeid = _m.elem.__uni10_typeid;
          elem.__ongpu = _m.elem.__ongpu;
          elem.__elemNum = _m.elem.__elemNum;
          elem.__elem = _m.elem.__elem;
          return *this;
        };

        uni10_type operator[](uni10_uint64 idx) const{
          uni10_error_msg(idx > this->elem.__elemNum, "%s", "The input index exceed the element number.");
          uni10_type val = elem.__elem[idx];
          return val;
        };

        /*********************  OPERATOR **************************/

        /// @brief Print out Block
        ///
        /// Prints out the elements of Block
        ///
        /// For a 2 x 3 matrix \c M,
        /// \code
        /// 2 x 3 = 6
        ///
        /// -0.254 -0.858 -0.447
        ///
        /// 0.392  0.331 -0.859
        /// \endcode
        /// The number of elements is 2 x 3 = 6 and followed by a 2 by 3 matrix of elements.
        friend std::ostream& operator<< <>(std::ostream& os, const Block& _b); // --> uni10_elem().print_elem()

         // Element-wise addition.
        template<typename _T>
          friend Matrix<_T> operator+(const Block<uni10_double64>& b1, const Block<_T>& b2);

        template<typename _T>
          friend Matrix<uni10_complex128> operator+(const Block<uni10_complex128>& b1, const Block<_T>& b2);

         // Element-wise subtract.
        template<typename _T>
          friend Matrix<_T> operator-(const Block<uni10_double64>& b1, const Block<_T>& b2);

        template<typename _T>
          friend Matrix<uni10_complex128> operator-(const Block<uni10_complex128>& b1, const Block<_T>& b2);

         // Element-wise multiplication.
        template<typename _T>
          friend Matrix<_T> operator*(const Block<uni10_double64>& b1, const Block<_T>& b2);

        template<typename _T>
          friend Matrix<uni10_complex128> operator*(const Block<uni10_complex128>& b1, const Block<_T>& b2);

        template<typename _T, typename _U>
          friend uni10_bool operator== (const Block<_T>& m1, const Block<_U>& m2);

        template<typename _T, typename _U>
          friend uni10_bool operator!= (const Block<_T>& m1, const Block<_U>& m2);

        //UNI10_LINALG_RETURN_VALUE
        template<typename _uni10_type> 
          friend Matrix<_uni10_type> getDiag( const Block<_uni10_type>& A );

        template<typename _To, typename _T, typename _U> 
          friend void dot( Matrix<_To>& C, const Block<_T>& A, const Block<_U>& B, UNI10_INPLACE on );

        template<typename _T> 
          friend void dot( Matrix<uni10_complex128>& A, const Block<_T>& B, UNI10_INPLACE on );

        friend void dot( Matrix<uni10_double64>& A, const Block<uni10_double64>& B, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > qr( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > rq( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > lq( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > ql( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > qdr( const Block<_uni10_type>& M);

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > ldq( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > qdr_cpivot( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > svd( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > sdd( const Block<_uni10_type>& M );

        template<typename _uni10_type>
          friend std::vector< Matrix<_uni10_type> > eigh( const Block<_uni10_type>& _Mij);

        template<typename _uni10_type>
          friend std::vector< Matrix<uni10_complex128> > eig( const Block<_uni10_type>& _Mij);

        template<typename _uni10_type>
          friend _uni10_type sum( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend uni10_double64 norm( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend Matrix<_uni10_type> inverse( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend Matrix<_uni10_type> transpose( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend Matrix<_uni10_type> dagger( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend Matrix<_uni10_type> conj( const Block<_uni10_type>& Mij );

        template<typename _uni10_type>
          friend _uni10_type det( const Block<_uni10_type>& _Mij );

        template<typename _uni10_type>
          friend _uni10_type trace( const Block<_uni10_type>& _Mij );

        template<typename _uni10_type>
          friend Matrix<_uni10_type> exph( uni10_double64 a, const Block<_uni10_type>& mat);

        template<typename _uni10_type>
          friend void svd( const Block<_uni10_type>& Mij, Matrix<_uni10_type>& U, Matrix<_uni10_type>& S, Matrix<_uni10_type>& VT, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void sdd( const Block<_uni10_type>& Mij, Matrix<_uni10_type>& U, Matrix<_uni10_type>& S, Matrix<_uni10_type>& VT, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend class Matrix;

        template<typename _uni10_type>
          friend class UniTensor;

    };

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, const Block<uni10_type>& _b){

      _b.elem.print_elem(_b.Rnum, _b.Cnum, _b.diag);

      os << "";

      return os;
    }

  template<typename T, typename U>
    uni10_bool operator== (const Block<T>& b1, const Block<U>& b2){

      if( (b1.Rnum != b2.Rnum) || (b1.Cnum != b2.Cnum) || (b1.diag != b2.diag) )
        return false;

      //Check real part 
     
      for(uni10_int i = 0; i < b1.elem.__elemNum; i++)
        if(UNI10_REAL(b1.elem.__elem[i] )- UNI10_REAL(b2.elem.__elem[i] )> 10E-12)
          return false;

      if(b1.elem.__uni10_typeid == 2)
        for(uni10_int i = 0; i < b1.elem.__elemNum; i++)
          if(UNI10_IMAG(b1.elem.__elem[i]) - UNI10_IMAG(b2.elem.__elem[i]) > 10E-12)
            return false;

      return true; 
    }

  template<typename T, typename U>
    uni10_bool operator!= (const Block<T>& b1, const Block<U>& b2){

      return !(b1 == b2); 

    }

};

#endif
