#ifndef __UNI10_LINALG_SCALAPACK_MPI_D_H__
#define __UNI10_LINALG_SCALAPACK_MPI_D_H__

#include "uni10/uni10_type.h"

namespace uni10{

  namespace uni10_linalg{

    // Blas 
    //
    //Done
    void vectorAdd(uni10_double64 a, uni10_double64* X, uni10_int32 incx, uni10_double64* Y, uni10_int32 incy, uni10_uint64 N);   // Y = a*X + Y

    //Done
    void vectorAdd(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y + X

    //Done
    void vectorSub(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y - X

    //Done
    void vectorMul(uni10_double64* Y, uni10_double64* X, uni10_uint64 N);   // Y = Y * X, element-wise multiplication;

    //Done
    void vectorScal(uni10_double64 a, uni10_double64* X, uni10_uint64 N);  // X = a * X

    void vectorExp(uni10_double64 a, uni10_double64* X, uni10_uint64 N);

    //Done
    uni10_double64 vectorSum(const uni10_double64* X, uni10_uint64 mpxnq, uni10_int inc);

    //Done
    uni10_double64 vectorNorm(const uni10_double64* A, const uni10_int* desc, uni10_uint64 M, uni10_uint64 N);

    //Done
    void matrixDot(const uni10_double64* A, const uni10_int* descA, const uni10_double64* B, const uni10_int* descB, 
        const uni10_int M, const uni10_int N, const uni10_int K, uni10_double64* C, uni10_int* descC);

    //Done
    void diagRowMul(uni10_double64* mat, uni10_double64* diag, uni10_int mp, uni10_int nq, uni10_int rh);

    //Done
    void diagColMul(uni10_double64* mat, uni10_double64* diag, uni10_int mp, uni10_int nq, uni10_int rh);

    //Done
    void setTranspose(const uni10_double64* A, const uni10_int* descA, uni10_int M, uni10_int N, uni10_double64* AT, uni10_int* descAT);

    //Done
    void setDagger(const uni10_double64* A, const uni10_int* descA, uni10_int M, uni10_int N, uni10_double64* AT, uni10_int* descAT);

    void setIdentity(uni10_double64* A, uni10_uint64 M, uni10_uint64 N);

    void trimatrixEigH(uni10_double64* D, uni10_double64* E, uni10_int32 N, 
        uni10_double64* z=NULL, uni10_int32 LDZ=1);

    // Lapack
    //
    /*Generate a set of row vectors which form a othonormal basis
     *For the incoming matrix "elem", the number of row <= the number of column, M <= N
     */
    //Done
    void matrixSVD(const uni10_double64* Mij_ori, const uni10_int blkrow, const uni10_int blkcol, const uni10_int* descMij, 
        uni10_int M, uni10_int N, uni10_double64* U, uni10_int* descU, uni10_double64* S, uni10_double64* vT, uni10_int* descvT);

    //Done
    void matrixSDD(uni10_double64* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_double64* U, uni10_double64* S, uni10_double64* vT);

    //Done
    void matrixQR(const uni10_double64* Mij_ori, const uni10_int blkrow, const uni10_int blkcol, const uni10_int* descMij,
        uni10_int M, uni10_int N, uni10_double64* Q, uni10_int* descQ, uni10_double64* R, uni10_int* descR);

    //Done
    void matrixRQ(const uni10_double64* Mij_ori, const uni10_int blkrow, const uni10_int blkcol, const uni10_int* descMij,
        uni10_int M, uni10_int N, uni10_double64* R, uni10_int* descR, uni10_double64* Q, uni10_int* descQ);

    //Done
    void matrixQL(const uni10_double64* Mij_ori, const uni10_int blkrow, const uni10_int blkcol, const uni10_int* descMij,
        uni10_int M, uni10_int N, uni10_double64* Q, uni10_int* descQ, uni10_double64* L, uni10_int* descL);

    //Done
    void matrixLQ(const uni10_double64* Mij_ori, const uni10_int blkrow, const uni10_int blkcol, const uni10_int* descMij,
        uni10_int M, uni10_int N, uni10_double64* L, uni10_int* descL, uni10_double64* Q, uni10_int* descQ);

    void matrixQDR(uni10_double64* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_double64* Q, uni10_double64* D, uni10_double64* R);

    void matrixLDQ(uni10_double64* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_double64* L, uni10_double64* D, uni10_double64* Q);

    void matrixQDRCPIVOT(uni10_double64* Mij_ori, uni10_int32 M, uni10_int32 N, uni10_double64* Q, uni10_double64* D, uni10_double64* R);

    void matrixInv(uni10_double64* A, uni10_int32 N);

    uni10_double64 matrixDet(uni10_double64* A, uni10_int32 N);

    //=====================================================================================//
    //
    void reshapeElem(uni10_double64* elem, uni10_uint64* transOffset);

    void eigSyDecompose(const uni10_double64* a, const uni10_int mp, const uni10_int nq, const uni10_int* desca, 
        uni10_int32 N, uni10_double64* w, uni10_double64* z, uni10_int* descz);

  };/* namespace uni10_linalg */

};/* namespace uni10 */

#endif
