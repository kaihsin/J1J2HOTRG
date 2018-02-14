/****************************************************************************
*  @file uni10_lapack_wrapper.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2016
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
*  @brief C wrapper functions for fortran BLAS and LAPACK libraries
*  @author Ying-Jer Kao, Yun-Hsuan Chou
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef __UNI10_BLACS_WRAPPER_MPI_H__
#define __UNI10_BLACS_WRAPPER_MPI_H__

#include <stdint.h>

#include <complex>

#include "uni10/uni10_type.h"
//BLACS C interface
extern "C"

{

  extern void   Cblacs_pinfo( uni10_int* mypnum, uni10_int* nprocs);

  extern void   Cblacs_get( uni10_int context, uni10_int request, uni10_int* value);

  extern uni10_int    Cblacs_gridinit( uni10_int* context, char * order, uni10_int np_row, uni10_int np_col);

  extern void   Cblacs_gridinfo( uni10_int context, uni10_int*  np_row, uni10_int* np_col, uni10_int*  my_row, uni10_int*  my_col);

  extern void   Cblacs_gridexit( uni10_int context);

  extern void   Cblacs_exit( uni10_int error_code);

  extern void   Cblacs_barrier( uni10_int icontxt, char* scope );

  extern void   Cdgsum2d(uni10_int icontxt, char* scope , char* top, uni10_int m, uni10_int n, double* A, uni10_int lda, uni10_int rdest, uni10_int cdest);
 
  //Sending
  extern void   Cdgesd2d(uni10_int icontxt, uni10_int m , uni10_int n, double* A, 
      uni10_int lda, uni10_int rdest, uni10_int cdest);

  extern void   Cdgebs2d(uni10_int icontxt, char * scope, char* top, uni10_int m, uni10_int n,
      double* A, uni10_int lda);

  extern void   Cdtrsd2d(uni10_int icontxt, char* uplo, char*diag, uni10_int m, uni10_int n,
      double* A, uni10_int lda, uni10_int rdest, uni10_int cdest);

  extern void   Cdtrbs2d(uni10_int icontxt, char* scope, char* top, char* uplo, char* diag,
      uni10_int m , uni10_int n, double* A, uni10_int lda);

  //Receving
  extern void   Cdgerv2d(uni10_int icontxt, uni10_int m, uni10_int n, double* A, 
      uni10_int lda, uni10_int rsrc, uni10_int csrc);

  extern void   Cdgebr2d(uni10_int icontxt, char* scope, char* top, uni10_int m, uni10_int n,
      double* A, uni10_int lda, uni10_int rsrc, uni10_int csrc);

  extern void   Cdtrrv2d(uni10_int icontxt, char* uplo, char* diag, uni10_int m, uni10_int n, double* A, uni10_int lda, uni10_int rsrc, uni10_int csrc);

  extern void   Cdtrbr2d(uni10_int icontxt, char* scope, char* top, char* uplo, char* diag,
      uni10_int m, uni10_int n, double* A, uni10_int lda, uni10_int rsrc, uni10_int csrc);

}

#endif
