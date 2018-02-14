/****************************************************************************
*  @file uni10_auxiliary.hpp
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2016
*    Natioanl Taiwan University
*    National Tsing-Hua University
*  
*    This file is designed for setting environment variables and allocating uni10_elems in general.
*
*  @endlicense
*  @author Yun-Hsuan Chou
*  @date 2014-05-06
*  @since
*
*****************************************************************************/

#ifndef __UNI10_AUXILIARY_H__
#define __UNI10_AUXILIARY_H__

#include <stdio.h>

#include "uni10/uni10_type.h"
#include "uni10/uni10_env_info.h"

namespace uni10{

  //extern env_type(env_info, _type) env_variables;

  void uni10_create(int argc=0, char** argv=NULL);

  void uni10_destroy();

  void uni10_print_env_info();

  bool uni10_func(load_uni10_rc, _type)( );

  //inline void get_real_memsize(uni10_unit64& real_memsize){

}

#endif
