#ifndef __UNI10_RESIZE_SCALAPACK_MPI_H__
#define __UNI10_RESIZE_SCALAPACK_MPI_H__

#include "uni10/uni10_type.h"
#include "uni10/uni10_scalapack_mpi/uni10_elem_scalapack_mpi.h"

namespace uni10{

  template <typename _uni10_type>
    void resize(uni10_elem_scalapack_mpi<_uni10_type>& Eout, uni10_uint64 _row, uni10_uint64 _col, 
        const uni10_elem_scalapack_mpi<_uni10_type>& Ein, const uni10_uint64 _Rnum, const uni10_uint64 _Cnum, uni10_const_bool& _fixHead);

  void resize_(uni10_elem_double64& Eout, const uni10_elem_double64& Ein, uni10_const_bool& _fixHead);

  void resize_(uni10_elem_complex128& Eout, const uni10_elem_complex128& Ein, uni10_const_bool& _fixHead);

  template <typename _uni10_type>
    void resize(uni10_elem_scalapack_mpi<_uni10_type>& Eout, uni10_uint64 _row, uni10_uint64 _col, 
        const uni10_elem_scalapack_mpi<_uni10_type>& Ein, const uni10_uint64 _Rnum, const uni10_uint64 _Cnum, uni10_const_bool& _fixHead){

      // Parameters, _row, col, _Rnum, _Cnum, are useless because Rnum and Cnum are saved in the class uni10_elem_scalapack_mpi.
      resize_(Eout, Ein, _fixHead);

    }



} // End of namespace.

#endif
