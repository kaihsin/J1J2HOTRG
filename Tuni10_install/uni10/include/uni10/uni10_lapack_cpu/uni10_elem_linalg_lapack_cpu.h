#ifndef __UNI10_ELEM_LINALG_LAPACK_CPU_H__
#define __UNI10_ELEM_LINALG_LAPACK_CPU_H__

#include "uni10/uni10_type.h"
#include "uni10/uni10_lapack_cpu/uni10_elem_lapack_cpu.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_tools_lapack_cpu.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_linalg_lapack_cpu_d.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_linalg_lapack_cpu_dz.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_linalg_lapack_cpu_z.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_resize_lapack_cpu.h"

typedef uni10::uni10_elem_lapack_cpu<uni10_double64>     uni10_elem_double64;
typedef uni10::uni10_elem_lapack_cpu<uni10_complex128>   uni10_elem_complex128;

namespace uni10{

  // Blas 
  //
  // UNI10_DOUBLE64
  void vectorAdd(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N);

  void vectorSub(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N);

  void vectorMul(uni10_elem_double64* Y, const uni10_elem_double64* X, const uni10_uint64* N);   // Y = Y * X, element-wise multiplication;

  void vectorScal(uni10_double64* a, uni10_elem_double64* X, uni10_uint64* N);   // X = a * X

  void vectorExp(uni10_double64* a, uni10_elem_double64* X, uni10_uint64* N);

  uni10_double64 vectorSum (const uni10_elem_double64* X, const uni10_uint64* N, uni10_int32* inc);

  uni10_double64 vectorNorm(const uni10_elem_double64* X, const uni10_uint64* N, uni10_int32* inc);

  void matrixAdd(uni10_elem_double64* A, uni10_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N );

  void matrixAdd(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C);

  void matrixSub(uni10_elem_double64* A, uni10_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N );

  void matrixSub(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C);

  void matrixMul(uni10_elem_double64* A, uni10_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N );

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* C);

  void matrixDot(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_double64* C);

  uni10_double64 matrixTrace(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N);

  void setTranspose(const uni10_elem_double64* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* AT);

  void setTranspose(uni10_elem_double64* A, uni10_uint64* M, uni10_uint64* N);

  void setDagger(const uni10_elem_double64* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_double64* AT);

  void setDagger(uni10_elem_double64* A, uni10_uint64* M, uni10_uint64* N);

  void setConjugate(const uni10_elem_double64* A, const uni10_uint64* N, uni10_elem_double64* A_conj);

  void setConjugate(uni10_elem_double64* A, uni10_uint64* N);

  void setDiag(uni10_elem_double64* _elem, const uni10_elem_double64* src_elem, const uni10_uint64* M, const uni10_uint64* N);

  void setIdentity(uni10_elem_double64* A, const uni10_bool* is_diag, const uni10_uint64* M, const uni10_uint64* N);

  void setNormalRand(uni10_elem_double64* A, const uni10_bool* is_diag, const uni10_uint64* M, const uni10_uint64* N, 
      const uni10_double64* mu, const uni10_double64* var, const uni10_int64* seed);

  void setUniformRand(uni10_elem_double64* A, const uni10_bool* is_diag, const uni10_uint64* M, const uni10_uint64* N, 
      const uni10_double64* up, const uni10_double64* dn, const uni10_int64* seed);

  void trimatrixEigh(uni10_elem_double64* D, uni10_elem_double64* E, uni10_uint64* N, 
      uni10_elem_double64* z=NULL, uni10_uint64* LDZ=NULL);

  // Blas 
  //
  // UNI10_COMPLEX128
  void vectorAdd(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N);

  void vectorSub(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N);

  void vectorMul(uni10_elem_complex128* Y, const uni10_elem_complex128* X, const uni10_uint64* N);   // Y = Y * X, element-wise multiplication;

  void vectorScal(uni10_complex128* a, uni10_elem_complex128* X, uni10_uint64* N);   // X = a * X

  void vectorExp(uni10_complex128* a, uni10_elem_complex128* X, uni10_uint64* N);

  uni10_complex128 vectorSum (const uni10_elem_complex128* X, const uni10_uint64* N, uni10_int32* inc);

  uni10_double64   vectorNorm(const uni10_elem_complex128* X, const uni10_uint64* N, uni10_int32* inc);

  void matrixAdd(uni10_elem_complex128* A, uni10_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N);

  void matrixAdd(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixSub(uni10_elem_complex128* A, uni10_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N);

  void matrixSub(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdiag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixMul(uni10_elem_complex128* A, uni10_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag,
      const uni10_uint64* M, const uni10_uint64* N );

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixDot(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C);

  uni10_complex128 matrixTrace(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N);

  void setTranspose(const uni10_elem_complex128* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* AT);

  void setTranspose(uni10_elem_complex128* A, uni10_uint64* M, uni10_uint64* N);

  void setDagger(const uni10_elem_complex128* A, const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* AT);

  void setDagger(uni10_elem_complex128* A, uni10_uint64* M, uni10_uint64* N);

  void setConjugate(const uni10_elem_complex128* A, const uni10_uint64* N, uni10_elem_complex128* A_conj);

  void setConjugate(uni10_elem_complex128* A, uni10_uint64* N);

  void setDiag(uni10_elem_complex128* _elem, const uni10_elem_complex128* src_elem, const uni10_uint64* M, const uni10_uint64* N);

  void setIdentity(uni10_elem_complex128* A, const uni10_bool* is_diag, const uni10_uint64* M, const uni10_uint64* N);

  void setNormalRand(uni10_elem_complex128* A, const uni10_bool* is_diag, const uni10_uint64* M, const uni10_uint64* N, 
      const uni10_double64* mu, const uni10_double64* var, const uni10_int64* seed);

  void setUniformRand(uni10_elem_complex128* A, const uni10_bool* is_diag, const uni10_uint64* M, const uni10_uint64* N, 
      const uni10_double64* up, const uni10_double64* dn, const uni10_int64* seed);

  void trimatrixEigh(uni10_elem_complex128* D, uni10_elem_complex128* E, uni10_uint64* N, 
      uni10_elem_complex128* z=NULL, uni10_uint64* LDZ=NULL);

  // Blas 
  //
  // MIX
  void vectorAdd(uni10_elem_complex128* Y, const uni10_elem_double64* X, const uni10_uint64* N);

  void vectorSub(uni10_elem_complex128* Y, const uni10_elem_double64* X, const uni10_uint64* N);

  void vectorMul(uni10_elem_complex128* Y, const uni10_elem_double64* X, const uni10_uint64* N);   // Y = Y * X, element-wise multiplication;

  void vectorScal(uni10_double64* a, uni10_elem_complex128* X, uni10_uint64* N);   // X = a * X

  void vectorExp(uni10_double64* a, uni10_elem_complex128* X, uni10_uint64* N);

  void matrixAdd(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixSub(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixMul(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixDot(const uni10_elem_double64* A, uni10_const_bool* Aisdag, const uni10_elem_complex128* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C);

  void matrixAdd(uni10_elem_complex128* A, uni10_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N );

  void matrixAdd(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixSub(uni10_elem_complex128* A, uni10_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N );

  void matrixSub(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  // Lack
  void matrixMul(uni10_elem_complex128* A, uni10_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N);

  void matrixMul(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, uni10_elem_complex128* C);

  void matrixDot(const uni10_elem_complex128* A, uni10_const_bool* Aisdag, const uni10_elem_double64* B, uni10_const_bool* Bisdag, 
      const uni10_uint64* M, const uni10_uint64* N, const uni10_uint64* K, uni10_elem_complex128* C);

  // LAPACK
  //
  //UNI10_DOUBLE64
  void matrixQR(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* R);

  void matrixRQ(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* R, uni10_elem_double64* Q);

  void matrixQL(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* L);

  void matrixLQ(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* L, uni10_elem_double64* Q);

  void matrixQDR(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* D, uni10_elem_double64* R);

  void matrixQDRCPIVOT(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* Q, uni10_elem_double64* D, uni10_elem_double64* R);

  void matrixLDQ(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* L, uni10_elem_double64* D, uni10_elem_double64* Q);

  void matrixSVD(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* U, uni10_elem_double64* S, uni10_elem_double64* vT);

  void matrixSDD(const uni10_elem_double64* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_double64* U, uni10_elem_double64* S, uni10_elem_double64* vT);

  void matrixEigh(const uni10_elem_double64* Mij_ori, uni10_const_bool *isMdiag, const uni10_uint64* N, 
      uni10_elem_double64* Eig, uni10_elem_double64* EigVec);

  void matrixEig(const uni10_elem_double64* Mij_ori, uni10_const_bool *isMdiag, const uni10_uint64* N, 
      uni10_elem_complex128* Eig, uni10_elem_complex128* EigVec);

  void matrixInv(const uni10_elem_double64* A, const uni10_uint64* N, uni10_const_bool* isdiag);

  uni10_double64 matrixDet(const uni10_elem_double64* A, const uni10_uint64* N, uni10_const_bool* isdiag);

  // LAPACK
  //
  //UNI10_COMPLEX128
  void matrixQR(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* R);

  void matrixRQ(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* R, uni10_elem_complex128* Q);

  void matrixQL(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* L);

  void matrixLQ(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* L, uni10_elem_complex128* Q);

  void matrixQDR(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* D, uni10_elem_complex128* R);

  void matrixQDRCPIVOT(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* Q, uni10_elem_complex128* D, uni10_elem_complex128* R);

  void matrixLDQ(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* L, uni10_elem_complex128* D, uni10_elem_complex128* Q);

  void matrixSVD(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* U, uni10_elem_complex128* S, uni10_elem_complex128* vT);

  void matrixSDD(const uni10_elem_complex128* Mij_ori, uni10_const_bool* isdiag, const uni10_uint64* M, const uni10_uint64* N, 
      uni10_elem_complex128* U, uni10_elem_complex128* S, uni10_elem_complex128* vT);

  void matrixEigh(const uni10_elem_complex128* Mij_ori, uni10_const_bool *isMdiag, const uni10_uint64* N, 
      uni10_elem_complex128* Eig, uni10_elem_complex128* EigVec);

  void matrixEig(const uni10_elem_complex128* Mij_ori, uni10_const_bool *isMdiag, const uni10_uint64* N, 
      uni10_elem_complex128* Eig, uni10_elem_complex128* EigVec);

  void matrixInv(const uni10_elem_complex128* A, const uni10_uint64* N, uni10_const_bool* isdiag);

  uni10_complex128 matrixDet(const uni10_elem_complex128* A, const uni10_uint64* N, uni10_const_bool* isdiag);

}

#endif

