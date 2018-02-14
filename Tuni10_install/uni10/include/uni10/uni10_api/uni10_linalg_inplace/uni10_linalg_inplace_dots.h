#ifndef __UNI10_LINALG_INPLACE_DOTS_H__
#define __UNI10_LINALG_INPLACE_DOTS_H__

#include <utility>

#include "uni10/uni10_api/Matrix.h"
#include "uni10_linalg_inplace_dot.h"

namespace uni10{

  // The driver fucntions of dots.
  inline void dot_dd(void* m3, const void* m1, const void* m2){

    dot( *((Matrix<double>*)m3), *((Matrix<double>*)m1), *((Matrix<double>*)m2), INPLACE);

  }

  inline void dot_dz(void* m3, const void* m1, const void* m2){

    dot( *((Matrix<uni10_complex128 >*)m3), *((Matrix<double>*)m1), *((Matrix<uni10_complex128 >*)m2), INPLACE);

  }

  inline void dot_zd(void* m3, const void* m1, const void* m2){

    dot( *((Matrix<uni10_complex128 >*)m3), *((Matrix<uni10_complex128 >*)m1), *((Matrix<double>*)m2), INPLACE);

  }

  inline void dot_zz(void* m3, const void* m1, const void* m2){

    dot( *((Matrix<uni10_complex128 >*)m3), *((Matrix<uni10_complex128 >*)m1), *((Matrix<uni10_complex128 >*)m2), INPLACE);

  }

  // Function pointers of vector of dots' driver functions.
  static void (*dot_driver[])(void* m3, const void* m1, const void* m2) = {dot_dd, dot_zd, dot_dz, dot_zz};

  // Dots with mix type.
  template<typename T>
    void dots_mix(Matrix<T>& mout, const std::vector< std::pair<void*, int> >& _mlist, UNI10_INPLACE on);

  // Matrices in the mlist have same type.
  template<typename T> 
    void dots_pure(Matrix<T>& mout, const std::vector< std::pair<void*, int> >& _mlist, UNI10_INPLACE on);


  // Dots with pure type (User API level).
  template<typename uni10_type> 
    void dots(Matrix<uni10_type>& mout, const std::vector< Matrix<uni10_type>* >& _mlist, UNI10_INPLACE on);

  template<typename T>
    void dots_mix(Matrix<T>& mout, const std::vector< std::pair<void*, int> >& _mlist, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace");
      uni10_error_msg(mout.typeID() == 1, "%s", 
          "There are complex matrices in the input arguments. Hence, the output matrix must be complex.");

      if(_mlist.size() == 1){

        if(_mlist[0].second == 1)
          mout = *(Matrix<uni10_double64>*)(_mlist[0].first);
        else
          mout = *(Matrix<uni10_complex128>*)(_mlist[0].first);
        return;

      }

      std::vector<std::pair<void*, int> > mlist = _mlist;

      if(mlist.size() == 2){

        int driver_type = (mlist[0].second - 1) + (2 * (mlist[1].second - 1));
        if(driver_type == 0){
          Matrix<uni10_double64> _mout;
          dot_driver[driver_type]((void*)&_mout, mlist[0].first, mlist[1].first);
          mout = _mout;
        }else
          dot_driver[driver_type]((void*)&mout, mlist[0].first, mlist[1].first);

        return;

      }

      std::vector<std::pair<void*, int> > sub_mlist;
      int driver_type;
      std::vector<Matrix<uni10_double64> > buf_r(mlist.size()/2);
      std::vector<Matrix<uni10_complex128> > buf_c(mlist.size()/2);

      int offset = (mlist.size() % 2 == 0) ? 0 : 1;  

      if(offset == 1)
        sub_mlist.push_back(mlist[0]);

      for(int b = 0 ; b < (int)mlist.size()/2; b++){
        driver_type = (mlist[2*b+offset].second -1) + (2 * (mlist[2*b+1+offset].second - 1));
        if(driver_type == 0){
          dot_driver[driver_type]((void*)&buf_r[b], mlist[2*b+offset].first, mlist[2*b+1+offset].first);
          sub_mlist.push_back(std::pair<void*, int>((void*)&buf_r[b], 1));
        }else{
          dot_driver[driver_type]((void*)&buf_c[b], mlist[2*b+offset].first, mlist[2*b+1+offset].first);
          sub_mlist.push_back(std::pair<void*, int>((void*)&buf_c[b], 2));
        }

      }

      dots_mix(mout, sub_mlist, on);

    }

  template<typename T> 
    void dots_pure(Matrix<T>& mout, const std::vector< std::pair<void*, int> >& _mlist, UNI10_INPLACE on){

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace");

      if(_mlist.size() == 1){

        mout = *(Matrix<T>*)(_mlist[0].first);
        return;

      }

      std::vector<std::pair<void*, int> > mlist = _mlist;
      int driver_type = mout.typeID() == 1 ? 0 : 3;

      if(mlist.size() == 2){

        dot_driver[driver_type]((void*)&mout, mlist[0].first, mlist[1].first);
        return;

      }

      std::vector<std::pair<void*, int> > sub_mlist;
      std::vector< Matrix<T> > buf(mlist.size()/2);

      int offset = (mlist.size() % 2 == 0) ? 0 : 1;  

      if(offset == 1)
        sub_mlist.push_back(mlist[0]);

      for(int b = 0 ; b < (int)mlist.size()/2; b++){
        dot_driver[driver_type]((void*)&buf[b], mlist[2*b+offset].first, mlist[2*b+1+offset].first);
        sub_mlist.push_back(std::pair<void*, int>((void*)&buf[b], mout.typeID()));
      }

      dots_pure(mout, sub_mlist, on);

    }

  template<typename uni10_type> 
    void dots(Matrix<uni10_type>& mout, const std::vector< Matrix<uni10_type>* >& _mlist, UNI10_INPLACE on) {

      uni10_error_msg(on != 1, "%s", "Setting a wrong flag of uni10_Inplace." );

      if(_mlist.size() == 1){

        mout = (*_mlist[0]);
        return;

      }

      std::vector< Matrix<uni10_type>* > mlist = _mlist;

      if(mlist.size() == 2){

          dot(mout, *mlist[0], *mlist[1], INPLACE);

        return;

      }

      std::vector< Matrix<uni10_type>* > sub_mlist;
      int driver_type;
      std::vector<Matrix<uni10_type> > buf_mat(mlist.size()/2);

      int offset = (mlist.size() % 2 == 0) ? 0 : 1;  

      if(offset == 1)
        sub_mlist.push_back(mlist[0]);

      for(int b = 0 ; b < (int)mlist.size()/2; b++){

        dot(buf_mat[b], *mlist[2*b+offset], *mlist[2*b+1+offset], INPLACE);
        sub_mlist.push_back(&buf_mat[b]);

      }

      dots(mout, sub_mlist, on);

    }

}

#endif
