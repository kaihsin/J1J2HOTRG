#include <math.h>

#include "uni10/uni10_error.h"
#include "uni10/uni10_api/Matrix.h"

namespace uni10{

  template <> template <>
    void Matrix<uni10_complex128>::init(const UELEM(uni10_elem, _package, _type)<uni10_double64>& _m){

      this->elem.init(this->Rnum, this->Cnum, this->diag, _m);

    };

  template <> template <>
    void Matrix<uni10_double64>::init(const UELEM(uni10_elem, _package, _type)<uni10_complex128>& _m){

      this->elem.init(this->Rnum, this->Cnum, this->diag,_m);

    };

  // Copy constructor of RotC or CtoR.
  template<> template<>
    Matrix<uni10_complex128>::Matrix(Matrix<uni10_double64> const& _m): Block<uni10_complex128>(_m.Rnum, _m.Cnum, _m.diag){

      this->init(_m.elem);

    }

  template<> template<>
    Matrix<uni10_double64>::Matrix(Matrix<uni10_complex128> const& _m): Block<uni10_double64>(_m.Rnum, _m.Cnum, _m.diag){

      this->init(_m.elem);

    }

  template<> template<>
    Matrix<uni10_complex128>::Matrix(Block<uni10_double64> const& _b): Block<uni10_complex128>(_b.Rnum, _b.Cnum, _b.diag){

      this->init(_b.elem);

    }

  template<> template<>
    Matrix<uni10_double64>::Matrix(Block<uni10_complex128> const& _b): Block<uni10_double64>(_b.Rnum, _b.Cnum, _b.diag){

      this->init(_b.elem);

    }

  //
  // Assignment operator
  template<> template<>
    Matrix<uni10_double64>& Matrix<uni10_double64>::operator=(const Matrix<uni10_complex128>& _m){

      this->Rnum = _m.Rnum;
      this->Cnum = _m.Cnum;
      this->diag = _m.diag;
      this->init(_m.elem);
      return *this;

    }

  template<> template<>
    Matrix<uni10_complex128>& Matrix<uni10_complex128>::operator=(const Matrix<uni10_double64>& _m){

      this->Rnum = _m.Rnum;
      this->Cnum = _m.Cnum;
      this->diag = _m.diag;
      this->init(_m.elem);
      return *this;

    }

  template<> template<>
    Matrix<uni10_double64>& Matrix<uni10_double64>::operator=(const Block<uni10_complex128>& _b){

      this->Rnum = _b.Rnum;
      this->Cnum = _b.Cnum;
      this->diag = _b.diag;
      this->init(_b.elem);
      return *this;

    }

  template<> template<>
    Matrix<uni10_complex128>& Matrix<uni10_complex128>::operator=(const Block<uni10_double64>& _b){

      this->Rnum = _b.Rnum;
      this->Cnum = _b.Cnum;
      this->diag = _b.diag;
      this->init(_b.elem);
      return *this;

    }

  Matrix<uni10_complex128>& operator+=(Matrix<uni10_complex128>& _m1, const Matrix<uni10_double64>& _m2){
    matrixAdd(&_m1.elem, &_m1.diag, &_m2.elem, &_m2.diag, &_m2.Rnum, &_m2.Cnum );
    _m1.diag = (_m1.diag && _m2.diag);
    return _m1;
  }
  
  Matrix<uni10_complex128>& operator-=(Matrix<uni10_complex128>& _m1, const Matrix<uni10_double64>& _m2){
    matrixSub(&_m1.elem, &_m1.diag, &_m2.elem, &_m2.diag, &_m2.Rnum, &_m2.Cnum );
    _m1.diag = (_m1.diag && _m2.diag);
    return _m1;
  }
  
  Matrix<uni10_double64>& operator-=(Matrix<uni10_double64>& _m1, Matrix<uni10_complex128> a){
    uni10_error_msg(true, "%s", "Can't use Matrix<double>::operator*=(T ) as [T = std::complex<double>] ");
    return _m1;
  }

  Matrix<uni10_complex128>& operator*=(Matrix<uni10_complex128>& _m1, const Matrix<uni10_double64>& _m2){            
    matrixMul(&_m1.elem, &_m1.diag, &_m2.elem, &_m2.diag, &_m2.Rnum, &_m2.Cnum );
    _m1.diag = (_m1.diag || _m2.diag);
    return _m1;
  }

};

