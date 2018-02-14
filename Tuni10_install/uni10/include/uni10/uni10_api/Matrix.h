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
#ifndef __UNI10_MATRIX_H__
#define __UNI10_MATRIX_H__

#include "uni10/uni10_api/Block.h"

namespace uni10{

  /// @class Matrix
  /// @brief The Matrix class defines a common matrix
  ///
  /// Matrix is an auxilliary class used to extract block information on UniTensor and perform linear algebra
  /// operations. A symmetric tensor contains Block's with corresponding Qnum's. Each block is a Matrix with
  /// tensor elements.
  ///
  /// The Matrix follows the C convention that the memory storage is row-major and indices start from 0.
  
  // Element-wise multiplication.
  //template<typename uni10_type>
  //  Matrix<uni10_type> operator*(const Block<uni10_type>& Ma, const Block<uni10_type>& Mb); 
  template<typename T>
    Matrix<T> operator*(uni10_double64 a, const Block<T>& b1); 

  template<typename T>
    Matrix<uni10_complex128> operator*(uni10_complex128 a, const Block<T>* b1);

  template<typename T>
    Matrix<T> operator*( const Block<T>& Mb, uni10_double64 a); 

  template<typename T>
    Matrix<uni10_complex128> operator*( const Block<T>& Mb, uni10_complex128 a); 

  // Element-wise addition.
  template<typename T>
    Matrix<T> operator+(const Block<uni10_double64>& b1, const Block<T>& b2);

  template<typename T>
    Matrix<uni10_complex128> operator+(const Block<uni10_complex128>& b1, const Block<T>& b2);

  // Element-wise subtract.
  template<typename T>
    Matrix<T> operator-(const Block<uni10_double64>& b1, const Block<T>& b2); 

  template<typename T>
    Matrix<uni10_complex128> operator-(const Block<uni10_complex128>& b1, const Block<T>& b2); 

  // Element-wise multipliation.
  template<typename T>
    Matrix<T> operator*( const Block<uni10_double64>& b1, const Block<T>& b2); 

  template<typename T>
    Matrix<uni10_complex128> operator*( const Block<uni10_complex128>& b1, const Block<T>& b2); 

  Matrix<uni10_complex128>& operator+=(Matrix<uni10_complex128>& _m1, const Matrix<uni10_double64>& _m2);

  Matrix<uni10_complex128>& operator-=(Matrix<uni10_complex128>& _m1, const Matrix<uni10_double64>& _m2);

  Matrix<uni10_complex128>& operator*=(Matrix<uni10_complex128>& _m1, uni10_complex128 a);

  Matrix<uni10_complex128>& operator*=(Matrix<uni10_complex128>& _m1, const Matrix<uni10_double64>& _m2);

  template <typename uni10_type>
    class Matrix:public Block<uni10_type> {

      private:

        void uni10_elem_free();

        void set_elem_null();

        void init(const uni10_type* elem = NULL);

        void init(const UELEM(uni10_elem, _package, _type)<uni10_type>& elem);     // pointer to a real matrix

        template<typename U>
          void init(const UELEM(uni10_elem, _package, _type)<U>& elem );

      public:
        ///
        ///@brief Default constructor
        ///
        explicit Matrix();

        /// @brief Create a Matrix
        ///
        /// Allocate memory of size <tt> Rnum * Cnum </tt> ( or <tt> min(Rnum, Cnum)</tt> if \c diag is \c true) for
        /// matrix elements and set the elements to zero
        /// @param _Rnum Number of Rows
        /// @param _Cnum Number of Columns
        /// @param _diag Set \c true for diagonal matrix, defaults to \c false
        explicit Matrix(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _diag=false);

        explicit Matrix(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_type* _src, bool _diag=false);

        explicit Matrix(const std::string& fname);

        /// @brief Create a Matrix
        ///
        /// @param Block<T> or Matrix<T>
        Matrix<uni10_type>(Block<uni10_type> const& _b);

        /// @overload
        template<typename U>
          Matrix<uni10_type>(Block<U> const& _b);
        /// @overload
        
        Matrix(Matrix const& _m);

        /// @overload
        template<typename U>
          Matrix<uni10_type>(Matrix<U> const& _b);

        ~Matrix();

        void assign(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag = false);

        void load(const std::string& fname);

        void setElem(const uni10_type* elem, bool src_ongpu = false);

        void setElem(const std::vector<uni10_type>& elem, bool src_ongpu = false);

        void setDiag(const uni10_bool _isdiag);

        void set_zeros();
        
        void identity();

        void randomize(char UorN='U', uni10_double64 dn_mu=-1, uni10_double64 up_var=1, uni10_int64 seed=uni10_clock);

        void orthoRand(char UorN='U', uni10_double64 dn_mu=-1, uni10_double64 up_var=1, uni10_int64 seed=uni10_clock);

        /// @brief Assigns to Matrix
        ///
        /// Assigns the content of \c mat to Matrix, replacing the original content by reallocating new memory
        /// fit for \c mat.
        /// @param _m Second Matrix
        Matrix& operator=(const Matrix& _m);

        template<typename U>
          Matrix& operator=(const Matrix<U>& _m);

        Matrix& operator=(const Block<uni10_type>& _b);

        template<typename U>
          Matrix& operator=(const Block<U>& _m);

        uni10_type& operator[](uni10_uint64 idx){
          return this->elem[idx];
        };

        Matrix<uni10_type>& operator+=(const Matrix<uni10_type>& _m);

        Matrix<uni10_type>& operator-=(const Matrix<uni10_type>& _m);

        /// @brief Multiply Matrix by a scalar and assign
        ///
        /// Performs matrix elem-wise multiplication of two matrices \c Ma and \c Mb. Store the results in a new matrix
        Matrix<uni10_type>& operator*=(const Matrix<uni10_type>& _m);

        /// @brief Multiply Matrix by a scalar and assign
        ///
        /// Performs element-wise multiplication with a Real scalar \c a .
        Matrix<uni10_type>& operator*=(uni10_double64 a);

        /// @brief Multiplication of two matrices
        ///
        /// Performs matrix elemwise multiplication of two matrices \c Ma and \c Mb. Store the results in a new matrix 
        //friend Matrix<uni10_type> operator* <>(const Block& Ma, const Block& Mb); // Elem-wise multiplication

        /// @brief  Multiplication of a matrix and a scalar
        ///
        /// Performs element-wise  multiplication of \c Ma with \c  Store the results in a new matrix
        template<typename _T>
          friend Matrix<_T> operator* (uni10_double64 a, const Block<_T>& Ma);

        /// @overload
        template<typename _T>
          friend Matrix<uni10_complex128> operator* (uni10_complex128 a, const Block<_T>& Ma);

        /// @overload
        template<typename _T>
          friend Matrix<_T> operator* (const Block<_T>& Ma, uni10_double64 a);

        /// @overload
        template<typename _T>
          friend Matrix<uni10_complex128> operator* (const Block<_T>& Ma, uni10_complex128 a);

        friend Matrix<uni10_complex128>& operator+=(Matrix<uni10_complex128>& _m1, const Matrix<uni10_double64>& _m2);

        friend Matrix<uni10_complex128>& operator-=(Matrix<uni10_complex128>& _m1, const Matrix<uni10_double64>& _m2);

        friend Matrix<uni10_complex128>& operator*=(Matrix<uni10_complex128>& _m1, uni10_complex128 a);

        friend Matrix<uni10_complex128>& operator*=(Matrix<uni10_complex128>& _m1, const Matrix<uni10_double64>& _m2);

        template<typename Res, typename Obj, typename... Args> 
          friend void dot_args(Res& _m1, const Obj& _m2, const Args&... args);

        template<typename _uni10_type> 
          friend void dots(Matrix<_uni10_type>& _m, const std::vector< Matrix<_uni10_type*> >& mats, UNI10_INPLACE on);

        template<typename _uni10_type> 
          friend void resize( Matrix<_uni10_type>& Mout , const Matrix<_uni10_type>& Min, uni10_uint64 row, uni10_uint64 col, UNI10_INPLACE on);

        template<typename _uni10_type>
          friend void qr( const Matrix<_uni10_type>& Mij, Matrix<_uni10_type>& Q, Matrix<_uni10_type>& R, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void rq( const Matrix<_uni10_type>& Mij, Matrix<_uni10_type>& R, Matrix<_uni10_type>& Q, UNI10_INPLACE on  );

        template<typename _uni10_type>
          friend void lq( const Matrix<_uni10_type>& Mij, Matrix<_uni10_type>& L, Matrix<_uni10_type>& Q, UNI10_INPLACE on  );

        template<typename _uni10_type>
          friend void ql( const Matrix<_uni10_type>& Mij, Matrix<_uni10_type>& L, Matrix<_uni10_type>& Q, UNI10_INPLACE on  );

        template<typename _uni10_type>
          friend void qdr( const Matrix<_uni10_type>& Mij, Matrix<_uni10_type>& Q, Matrix<_uni10_type>& D, Matrix<_uni10_type>& R, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void ldq( const Matrix<_uni10_type>& Mij, Matrix<_uni10_type>& L, Matrix<_uni10_type>& D, Matrix<_uni10_type>& Q, UNI10_INPLACE on  );

        template<typename _uni10_type>
          friend void qdr_cpivot( const Matrix<_uni10_type>& Mij, Matrix<_uni10_type>& Q, Matrix<_uni10_type>& D, Matrix<_uni10_type>& R, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void eigh( const Matrix<_uni10_type>& Mij, Matrix<uni10_double64>& Eig, Matrix<_uni10_type>& EigVec, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void eig( const Matrix<_uni10_type>& Mij, Matrix<uni10_complex128>& Eig, Matrix<uni10_complex128>& EigVec, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void inverse( Matrix<_uni10_type>& Mij, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void transpose( Matrix<_uni10_type>& Mij, const Matrix<_uni10_type>& ori_Mij, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void transpose( Matrix<_uni10_type>& Mij, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void dagger( Matrix<_uni10_type>& Mij, const Matrix<_uni10_type>& ori_Mij, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void dagger( Matrix<_uni10_type>& Mij, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void conj( Matrix<_uni10_type>& Mij, const Matrix<_uni10_type>& ori_Mij, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void conj( Matrix<_uni10_type>& Mij, UNI10_INPLACE on );

    };

  template<typename T>
    Matrix<T> operator* (uni10_double64 a, const Block<T>& m1){
      Matrix<T> m2(m1);
      vectorScal(&a, &m2.elem, &m2.elem.__elemNum);
      return m2;
    }

  template<typename T>
    Matrix<uni10_complex128> operator* (uni10_complex128 a, const Block<T>& m1){
      Matrix<uni10_complex128> m2(m1);
      vectorScal(&a, &m2.elem, &m2.elem.__elemNum);
      return m2;
    }

  template<typename T>
    Matrix<T> operator* (const Block<T>& m1, uni10_double64 a){
      return a * m1;
    }

  template<typename T>
    Matrix<uni10_complex128> operator* (const Block<T>& m1, uni10_complex128 a){
      return a * m1;
    }

  template<typename T>
    Matrix<T> operator+ (const Block<uni10_double64>& m1, const Block<T>& m2){

      uni10_error_msg(m1.Rnum != m2.Rnum || m1.Cnum != m2.Cnum, "%s", "Lack err msg!!!");

      Matrix<T> m3(m1.Rnum, m1.Cnum, m1.diag && m2.diag);
      matrixAdd(&m1.elem, &m1.diag, &m2.elem, &m2.diag, &m1.Rnum, &m1.Cnum, &m3.elem);

      return m3;

    }
 
  template<typename T>
    Matrix<uni10_complex128> operator+ (const Block<uni10_complex128>& m1, const Block<T>& m2){

      uni10_error_msg(m1.Rnum != m2.Rnum || m1.Cnum != m2.Cnum, "%s", "Lack err msg!!!");

      Matrix<uni10_complex128> m3(m1.Rnum, m1.Cnum, m1.diag && m2.diag);
      matrixAdd(&m1.elem, &m1.diag, &m2.elem, &m2.diag, &m1.Rnum, &m1.Cnum, &m3.elem);

      return m3;

    }

  template<typename T>
    Matrix<T> operator- (const Block<uni10_double64>& m1, const Block<T>& m2){

      uni10_error_msg(m1.Rnum != m2.Rnum || m1.Cnum != m2.Cnum, "%s", "Lack err msg!!!");

      Matrix<T> m3(m1.Rnum, m1.Cnum, m1.diag && m2.diag);
      matrixSub(&m1.elem, &m1.diag, &m2.elem, &m2.diag, &m1.Rnum, &m1.Cnum, &m3.elem);

      return m3;
       
    }

  template<typename T>
    Matrix<uni10_complex128> operator- (const Block<uni10_complex128>& m1, const Block<T>& m2){

      uni10_error_msg(m1.Rnum != m2.Rnum || m1.Cnum != m2.Cnum, "%s", "Lack err msg!!!");

      Matrix<uni10_complex128> m3(m1.Rnum, m1.Cnum, m1.diag && m2.diag);
      matrixSub(&m1.elem, &m1.diag, &m2.elem, &m2.diag, &m1.Rnum, &m1.Cnum, &m3.elem);

      return m3;
       
    }

  template<typename T>
    Matrix<T> operator* (const Block<uni10_double64>& m1, const Block<T>& m2){

      uni10_error_msg(m1.Rnum != m2.Rnum || m1.Cnum != m2.Cnum, "%s", "Lack err msg!!!");

      Matrix<T> m3(m1.Rnum, m1.Cnum, m1.diag || m2.diag);
      matrixMul(&m1.elem, &m1.diag, &m2.elem, &m2.diag, &m1.Rnum, &m1.Cnum, &m3.elem);

      return m3;
       
    }

  template<typename T>
    Matrix<uni10_complex128> operator* (const Block<uni10_complex128>& m1, const Block<T>& m2){

      uni10_error_msg(m1.Rnum != m2.Rnum || m1.Cnum != m2.Cnum, "%s", "Lack err msg!!!");

      Matrix<uni10_complex128> m3(m1.Rnum, m1.Cnum, m1.diag || m2.diag);
      matrixMul(&m1.elem, &m1.diag, &m2.elem, &m2.diag, &m1.Rnum, &m1.Cnum, &m3.elem);

      return m3;
       
    }


};  /* namespace uni10 */

#endif /* MATRIX_H */
