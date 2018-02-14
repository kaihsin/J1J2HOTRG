#ifndef __UNI10_LINALG_SCALAPACK_MPI_Z_H__
#define __UNI10_LINALG_SCALAPACK_MPI_Z_H__

#include "uni10/uni10_type.h"

namespace uni10{

  namespace uni10_linalg{

    // Blas
    //
    void vectorAdd(uni10_complex128 a, uni10_complex128* X, uni10_int32 incx, uni10_complex128* Y, uni10_int32 incy, uni10_uint64 N);   // Y = a*X + Y

    void vectorAdd(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N);// Y = Y + X

    void vectorSub(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N);// Y = Y - X

    void vectorMul(uni10_complex128* Y, uni10_complex128* X, uni10_uint64 N); // Y = Y * X, element-wise multiplication;

    void vectorScal(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);// X = a * X

    void vectorExp(uni10_complex128 a, uni10_complex128* X, uni10_uint64 N);

    uni10_complex128 vectorSum(uni10_complex128* X, uni10_uint64 N, uni10_int32 inc);

    uni10_double64 vectorNorm(uni10_complex128* X, uni10_uint64 N, uni10_int32 inc);

    void matrixDot(uni10_complex128* A, uni10_complex128* B, uni10_int32 M, uni10_int32 N, uni10_int32 K, uni10_complex128* C);

    void diagRowMul(uni10_complex128* mat, uni10_complex128* diag, uni10_uint64 M, uni10_uint64 N);

    void diagColMul(uni10_complex128* mat, uni10_complex128* diag, uni10_uint64 M, uni10_uint64 N);

    void setTranspose(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N, uni10_complex128* AT);

    void setTranspose(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    void setDagger(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N, uni10_complex128* AT);

    void setDagger(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    void setConjugate(uni10_complex128 *A, uni10_uint64 N, uni10_complex128 *A_conj);

    void setConjugate(uni10_complex128 *A, uni10_uint64 N);

    void setIdentity(uni10_complex128* A, uni10_uint64 M, uni10_uint64 N);

    // Lapack
    //
    void matrixSVD(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* U, uni10_complex128* S, uni10_complex128* vT);

    void matrixSDD(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* U, uni10_complex128* S, uni10_complex128* vT);

    void matrixQR(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* Q, uni10_complex128* R);

    void matrixRQ(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* R, uni10_complex128* Q);

    void matrixQL(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* Q, uni10_complex128* L);

    void matrixLQ(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* L, uni10_complex128* Q);

    void matrixQDR(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* Q, uni10_complex128* D, uni10_complex128* R);

    void matrixLDQ(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* L, uni10_complex128* D, uni10_complex128* Q);

    void matrixQDRCPIVOT(uni10_complex128* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_complex128* Q, uni10_complex128* D, uni10_complex128* R);

    void matrixInv(uni10_complex128* A, uni10_int32 N);

    uni10_complex128 matrixDet(uni10_complex128* A, uni10_int32 N);

    //=================================================================================//

    void eigDecompose(uni10_complex128* Kij, uni10_int32 N, uni10_complex128* Eig, uni10_complex128 *EigVec);


  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
