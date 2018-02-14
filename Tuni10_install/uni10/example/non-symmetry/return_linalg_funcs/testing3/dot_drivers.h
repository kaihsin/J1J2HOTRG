#ifndef __UNI10_DOT_DRIVERS_H__
#define __UNI10_DOT_DRIVERS_H__

#include "uni10/uni10.hpp"

using namespace std;
using namespace uni10;

inline void dot_dd(void* m3, const void* m1, const void* m2){

  dot( *((Matrix<double>*)m3), *((Matrix<double>*)m1), *((Matrix<double>*)m2), INPLACE);

}

inline void dot_dz(void* m3, const void* m1, const void* m2){

  dot( *((Matrix<complex<double> >*)m3), *((Matrix<double>*)m1), *((Matrix<complex<double> >*)m2), INPLACE);

}

inline void dot_zd(void* m3, const void* m1, const void* m2){

  dot( *((Matrix<complex<double> >*)m3), *((Matrix<complex<double> >*)m1), *((Matrix<double>*)m2), INPLACE);

}

inline void dot_zz(void* m3, const void* m1, const void* m2){

  dot( *((Matrix<complex<double> >*)m3), *((Matrix<complex<double> >*)m1), *((Matrix<complex<double> >*)m2), INPLACE);

}


static void (*dot_driver[])(void* m3, const void* m1, const void* m2) = {dot_dd, dot_zd, dot_dz, dot_zz};


#endif
