#ifndef __UNI10_ELEM_SCALAPACK_MPI_H__
#define __UNI10_ELEM_SCALAPACK_MPI_H__

#include <iostream>

#include "uni10/uni10_env_info/uni10_env_info_scalapack_mpi.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_tools_scalapack_mpi.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_scalapack_tools_wrapper_mpi.h"

namespace uni10{

  template<typename uni10_type>
    class uni10_elem_scalapack_mpi{
      
      public:

        // Done
        explicit uni10_elem_scalapack_mpi();

        // Done
        explicit uni10_elem_scalapack_mpi(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag = false, uni10_bool _dist = false);

        // Done
        explicit uni10_elem_scalapack_mpi(const uni10_type* src, uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag = false, uni10_bool _dist = false);

        // Done
        explicit uni10_elem_scalapack_mpi(const uni10_elem_scalapack_mpi& _elem);

        // Done
        uni10_elem_scalapack_mpi& operator=(const uni10_elem_scalapack_mpi& _m){
          __uni10_typeid = _m.__uni10_typeid;
          __ongpu        = _m.__ongpu;
          this->copy(_m);
          return *this;
        }

        uni10_type& operator[](const uni10_uint64 idx){
          uni10_error_msg(idx>this->__elemNum, "%s", "The index is exceed the number of elements");
          return this->__elem[idx];
        }

        // Done
        ~uni10_elem_scalapack_mpi();

        // Done
        inline bool empty() const{ return status; }

        // Done
        void set_zeros();

        // Done
        void setElem(const uni10_type* src, bool src_ongpu = false);

        // Done
        void assign(uni10_uint64& _Rnum, uni10_uint64& _Cnum);

        // Done
        void clear();

        // Done
        void init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, const uni10_type* src=NULL);

        // Done
        void copy(const uni10_elem_scalapack_mpi& _elem);

        void copy(uni10_uint64 begin_idx, const uni10_elem_scalapack_mpi<uni10_type>& src, uni10_uint64 begin_src_idx ,uni10_uint64 len);

        void copy(uni10_uint64 begin_idx, const uni10_elem_scalapack_mpi<uni10_type>& src, uni10_uint64 len);

        void catElem(const uni10_elem_scalapack_mpi<uni10_type>& src);

        void print_elem(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag) const;

        void save(FILE* fp) const;

        void load(FILE* fp);

        uni10_type_id __uni10_typeid;  

        uni10_bool __ongpu;

        uni10_uint64 __elemNum;

        uni10_type* __elem;

        //
        // Additional parameters for controling memory in mpi version.
        //
        uni10_int ictxt, rank, npnum, nprow, npcol, blockrow, blockcol, rgrid, cgrid, myrow, mycol;

        //
        // Please put C elements into a tranposed matrix descriptor because scalapack is read by column ( Fortran )
        //
        uni10_int Rnum_, Cnum_;

        //
        // The heads of each blocks
        // Notice that the dists is counted by row (C/C++)
        // Only allocate on the master process.
        uni10_int r_head, r_offset, c_head, c_offset;
        uni10_int *r_lens, *c_lens;
        uni10_int *rdists, *cdists, *gdists;

        //
        // Descriptor for BLACS.
        //
        uni10_int desc[9];

        //
        // Master process.
        //
        uni10_int master;

        //
        // Is initialization or not.
        //
        uni10_int status;

      private:

        void init_dists();

    };

}

#endif
