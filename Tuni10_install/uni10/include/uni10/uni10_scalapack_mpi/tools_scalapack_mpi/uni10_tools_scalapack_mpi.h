#ifndef __UNI10_TOOLS_SCALAPACK_MPI_H__
#define __UNI10_TOOLS_SCALAPACK_MPI_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mpi.h"

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"

namespace uni10{

  void* uni10_elem_alloc(uni10_uint64 memsize);

  void* uni10_elem_copy(void* des, const void* src, uni10_uint64 memsize);

  void uni10_elem_free(void* ptr, uni10_uint64 memsize);

  void uni10_elemBzero(void* ptr, uni10_uint64 memsize);

  //For double ptr.
  //
  void uni10_setDiag(uni10_double64* elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_getDiag(uni10_double64* elem, uni10_double64* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_print_elem_i(const uni10_double64& elem_i);

  void uni10_getUpTri(uni10_double64* elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n); // elem -> tri_elem

  void uni10_getDnTri(uni10_double64* elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n);
  //For complex ptr.
  //
  void uni10_setDiag(uni10_complex128* elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_getDiag(uni10_complex128* elem, uni10_complex128* diag_elem, uni10_uint64 M, uni10_uint64 N, uni10_uint64 diag_N);

  void uni10_print_elem_i(const uni10_complex128& elem_i);
  
  void uni10_getUpTri(uni10_complex128* elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n); // elem -> tri_elem

  void uni10_getDnTri(uni10_complex128* elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n);
  // Convert
  void uni10_elem_cast(uni10_complex128* des, uni10_double64* src, uni10_uint64 N);

  void uni10_elem_cast(uni10_double64 *des, uni10_complex128 *src, uni10_uint64 N);

  void ToReal(uni10_double64& M_i, uni10_double64 val);

  void ToReal(uni10_complex128& M_i, uni10_double64 val);

  void ToComplex(uni10_double64& M_i, uni10_double64 val);

  void ToComplex(uni10_complex128& M_i, uni10_double64 val);

  void shrinkWithoutFree(uni10_uint64 memsize);

  void printDesc(const uni10_int* desc);

  void broadcast(uni10_double64* dest, const uni10_double64* buffer, const uni10_int count , const uni10_int root);

  void broadcast(uni10_complex128* dest, const uni10_complex128* buffer, const uni10_int count , const uni10_int root);

  void split_mast2dist(const uni10_double64* src, const uni10_int* r_lens, const uni10_int& nprow, const uni10_int* c_lens, const uni10_int& npcol, 
      const uni10_int* gdists, const uni10_int& Rnum_, const uni10_int& src_root, 
      uni10_double64* dist, const uni10_int& blkrow, const uni10_int& blkcol, const uni10_int& rank);

  void split_mast2dist(const uni10_complex128* src, const uni10_int* r_lens, const uni10_int& nprow, const uni10_int* c_lens, const uni10_int& npcol, 
      const uni10_int* gdists, const uni10_int& Rnum_, const uni10_int& src_root, 
      uni10_complex128* dist, const uni10_int& blkrow, const uni10_int& blkcol, const uni10_int& rank);

  void gather_dist2mast(const uni10_double64* dist, const uni10_int& blkrow, const uni10_int& blkcol, const uni10_int& rank, 
      uni10_double64* dest_ptr, const uni10_int* r_lens, const uni10_int& nprow, const uni10_int* const c_lens, const uni10_int& npcol, 
      const uni10_int* gdists, const uni10_int& Rnum_, const uni10_int& dest);

  void gather_dist2mast(const uni10_complex128* dist, const uni10_int& blkrow, const uni10_int& blkcol, const uni10_int& rank, 
      uni10_complex128* dest_ptr, const uni10_int* r_lens, const uni10_int& nprow, const uni10_int* c_lens, const uni10_int& npcol, 
      const uni10_int* gdists, const uni10_int& Rnum_, const uni10_int& dest);

}/* namespace uni10 */
#endif /* UNI10_TOOLS_H */
