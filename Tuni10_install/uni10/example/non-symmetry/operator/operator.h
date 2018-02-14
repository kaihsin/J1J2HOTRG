#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include "uni10/uni10.hpp"

void spin_check(uni10_float32 spin);

uni10::Matrix<uni10_double64> matSx(uni10_float32 spin=0.5);

uni10::Matrix<uni10_double64> matSy(uni10_float32 spin=0.5);

uni10::Matrix<uni10_double64> matSp(uni10_float32 spin=0.5);

uni10::Matrix<uni10_double64> matSm(uni10_float32 spin=0.5);

uni10::Matrix<uni10_double64> matSz(uni10_float32 spin=0.5);

uni10::UniTensor<uni10_double64> periodicSx(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

uni10::UniTensor<uni10_double64> periodicSy(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

uni10::UniTensor<uni10_double64> periodicSz(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

uni10::UniTensor<uni10_double64> periodicsqSx(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

uni10::UniTensor<uni10_double64> periodicsqSy(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

uni10::UniTensor<uni10_double64> periodicsqSz(uni10_uint64 siteNum, uni10_float32 spin = 0.5);

#endif
