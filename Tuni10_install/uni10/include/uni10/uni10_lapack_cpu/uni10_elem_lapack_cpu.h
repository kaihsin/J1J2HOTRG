#ifndef __UNI10_ELEM_LAPACK_CPU_H__
#define __UNI10_ELEM_LAPACK_CPU_H__

#include <iostream>
#include <vector>

#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_lapack_cpu/tools_lapack_cpu/uni10_tools_lapack_cpu.h"

namespace uni10{

  template<typename uni10_type>
    class uni10_elem_lapack_cpu{
      
      public:

        explicit uni10_elem_lapack_cpu();

        explicit uni10_elem_lapack_cpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag = false, uni10_bool _ongpu = false);

        uni10_elem_lapack_cpu(uni10_elem_lapack_cpu const& _elem);

        template<typename U>
          uni10_elem_lapack_cpu<uni10_type>(uni10_elem_lapack_cpu<U> const& _elem);

        uni10_elem_lapack_cpu& operator=(uni10_elem_lapack_cpu const& _m);

        template<typename U>
          uni10_elem_lapack_cpu<uni10_type>& operator=(uni10_elem_lapack_cpu<U> const& _elem);

        uni10_type& operator[](const uni10_uint64 idx){
          uni10_error_msg(idx>this->__elemNum, "%s", "The index is exceed the number of elements");
          return this->__elem[idx];
        }

        ~uni10_elem_lapack_cpu();

        inline bool empty() const{ return __elem == NULL; }

        void set_zeros();

        void setElem(const uni10_type* src, bool src_ongpu = false);

        void clear();

        void init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, uni10_elem_lapack_cpu const& _m);

        template<typename U>
          void init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, uni10_elem_lapack_cpu<U> const& _m);

        void init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, const uni10_type* src=NULL);

        void copy(uni10_elem_lapack_cpu const& _elem);

        template<typename U>
          void copy(uni10_elem_lapack_cpu<U> const& _elem);

        void copy(uni10_uint64 begin_idx, const uni10_elem_lapack_cpu<uni10_type>& src, uni10_uint64 begin_src_idx ,uni10_uint64 len);

        void copy(uni10_uint64 begin_idx, const uni10_elem_lapack_cpu<uni10_type>& src, uni10_uint64 len);

        void catElem(const uni10_elem_lapack_cpu<uni10_type>& src);

        void print_elem(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag) const;

        void save(FILE* fp) const;

        void load(FILE* fp);

        void reshape(const std::vector<int>& ori_bdDims, const std::vector<int>& new_bdIdx, uni10_elem_lapack_cpu& out_elem, bool inorder);

        uni10_type_id __uni10_typeid;  

        uni10_bool __ongpu;

        uni10_uint64 __elemNum;

        uni10_type* __elem;

    };

}

#endif
