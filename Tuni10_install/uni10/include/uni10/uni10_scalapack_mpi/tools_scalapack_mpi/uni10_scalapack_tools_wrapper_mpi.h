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
#ifndef __UNI10_SCALAPACK_TOOLS_WRAPPER_MPI_H__
#define __UNI10_SCALAPACK_TOOLS_WRAPPER_MPI_H__

#include <stdint.h>

#include <complex>

#include "uni10/uni10_type.h"
//BLACS C interface
extern "C"

{

  uni10_int numroc_( uni10_int *n, uni10_int *nb, uni10_int *iproc, uni10_int *isrcproc, uni10_int *nprocs);

  uni10_int indxg2p_( uni10_int* indxglob, uni10_int* nb, uni10_int* iproc, uni10_int* isrcproc, uni10_int* nprocs);

  void descinit_( uni10_int *desc, uni10_int *m, uni10_int *n, uni10_int *mb, 
      uni10_int *nb, uni10_int *irsrc, uni10_int *icsrc,
      uni10_int *ictxt, uni10_int *lld, uni10_int *info);

}

inline uni10_int numroc( uni10_int *n, uni10_int *nb, uni10_int *iproc, uni10_int *isrcproc, uni10_int *nprocs){
  return numroc_( n, nb, iproc, isrcproc, nprocs);
}

inline uni10_int indxg2p( uni10_int* indxglob, uni10_int* nb, uni10_int* iproc, uni10_int* isrcproc, uni10_int* nprocs){
  return indxg2p_(indxglob, nb, iproc, isrcproc, nprocs);
}

inline void descinit( uni10_int *desc, uni10_int *m, uni10_int *n, uni10_int *mb, 
      uni10_int *nb, uni10_int *irsrc, uni10_int *icsrc,
      uni10_int *ictxt, uni10_int *lld, uni10_int *info){
  descinit_(desc, m, n, mb, nb, irsrc, icsrc, ictxt, lld, info);
}

#endif
