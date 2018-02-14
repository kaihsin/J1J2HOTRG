#ifndef __UNI10_TYPE_H__
#define __UNI10_TYPE_H__

#include <stdio.h>
#include <string.h>

#include <complex>

namespace uni10{

  enum exu_type{
    no    =   0,
    cpu   =   1,
    gpu   =   2,
    mpi   =   3,
    mgpu  =   4
  };

};

// The excute type of Uni10.
#ifdef UNI_CPU
#define _type _cpu
#elif  UNI_GPU
#define _type _gpu
#elif  UNI_MPI
#define _type _mpi
#endif

#ifdef UNI_LAPACK
#define   _package   _lapack
#elif  UNI_CUSOLVER
#define   _package   _cusolver
#elif  UNI_SCALAPACK
#define   _package   _scalapack
#endif


#define uni10_clock            -1

#if defined(UNI_MKL)
#include "mkl_types.h"
typedef MKL_INT                uni10_int;
#define MKL_Complex8  std::complex<float>
#define MKL_Complex16 std::complex<double>
#else
typedef int                    uni10_int;
#endif

typedef int                    uni10_int32;       // shout int.
typedef long long int          uni10_int64;       // MKL uses long long int, not int64_t.
typedef size_t                 uni10_uint64;      
typedef float                  uni10_float32;
typedef double                 uni10_double64;
typedef std::complex<float>    uni10_complex64;
typedef std::complex<double>   uni10_complex128;

typedef int                    uni10_exu_type;    // To store the exu_type.
typedef int                    uni10_type_id;     // To store the typeid of the objects.

typedef bool                   uni10_bool;   
typedef const bool             uni10_const_bool;  

// Generate the typename of the uni10 system information.
#define env_type_helper(envinfo, type)  envinfo##type
#define env_type(envinfo, type)  env_type_helper(envinfo, type)

// Generate the suitable function name.
#define uni10_func_helper(func, type)  func##type
#define uni10_func(func, type)  uni10_func_helper(func, type)

// Generate the class name of uni10_elem.
#define elem_type_helper(UNI10ELEM, type)  UNI10ELEM##type
#define elem_type(UNI10ELEM, type)  elem_type_helper(UNI10ELEM, type)

#define UELEM_helper(UNI10ELEM, package, type) UNI10ELEM##package##type
#define UELEM(UNI10ELEM, package, type) UELEM_helper(UNI10ELEM, package, type)
//typedef double*               uni10_double_ptr;
//typedef std::complex<double>* uni10_complex_ptr;
//
//
#define UNI10_TYPE_ID(x)     (sizeof(x)/8)

#define Z_REAL(x)       (x).real()
#define Z_IMAG(x)       (x).imag()

inline uni10_double64 UNI10_REAL(uni10_double64 a){

  return a;

}

inline uni10_double64 UNI10_IMAG(uni10_double64 a){

  return 0.;

}
inline uni10_double64 UNI10_REAL(uni10_complex128 a){

  return a.real();

}

inline uni10_double64 UNI10_IMAG(uni10_complex128 a){

  return a.imag();

}

#endif
