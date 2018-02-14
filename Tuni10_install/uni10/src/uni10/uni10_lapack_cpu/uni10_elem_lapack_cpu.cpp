#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_lapack_cpu/uni10_elem_lapack_cpu.h"

namespace uni10{

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(): __uni10_typeid(UNI10_TYPE_ID(uni10_type)), __ongpu(false), __elemNum(0), __elem(NULL){};

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, uni10_bool _ongpu): 
      __uni10_typeid(UNI10_TYPE_ID(uni10_type)), __ongpu(_ongpu),__elem(NULL){

      init(_Rnum, _Cnum, _isdiag, NULL);

    }

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::uni10_elem_lapack_cpu(uni10_elem_lapack_cpu<uni10_type> const& _elem): 
      __uni10_typeid(_elem.__uni10_typeid), __ongpu(_elem.__ongpu), __elemNum(_elem.__elemNum), __elem(NULL){

      uni10_error_msg(_elem.__elem != NULL && _elem.__elemNum==0, "%s", "The pointer of this element container point to another pointer and does not allocate. Hence we can't copy it.");
      init(1, _elem.__elemNum, false, _elem);

    };

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>::~uni10_elem_lapack_cpu(){

      if(__elem != NULL && __elemNum != 0)
        uni10_elem_free(__elem, __elemNum * sizeof(uni10_type));

      __elem    = NULL;
      __elemNum = 0;

    };

  template<typename uni10_type>
    uni10_elem_lapack_cpu<uni10_type>& uni10_elem_lapack_cpu<uni10_type>::operator=(uni10_elem_lapack_cpu const& _m){

      __uni10_typeid = _m.__uni10_typeid;
      __ongpu        = _m.__ongpu;
      this->init(1, _m.__elemNum, false, _m.__elem);
      return *this;

    }

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::set_zeros(){

      uni10_error_msg( __elem == NULL || __elemNum == 0, "%s", "Please initialize the uni10_elem with the constructor uni10(uni10_uint64, uni10_uint64, bool) befero setting the elements.");

      uni10_uint64 memsize = __elemNum * sizeof(uni10_type);

      uni10_elemBzero( __elem, memsize );

    };

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::setElem(const uni10_type* src, bool src_ongpu){

      uni10_error_msg( src_ongpu, "%s", " The source pointer is on the device. Please install MAGMA or CUDAONLY gpu version instead.");
      uni10_error_msg( __elem == NULL, "%s", "Please initialize the uni10_elem with the constructor uni10(uni10_uint64, uni10_uint64, bool) befero setting the elements.");
      uni10_error_msg( src  == NULL, "%s", "The source ptr is NULL.");

      uni10_elem_copy( __elem, src, __elemNum * sizeof(uni10_type) );

    };

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, uni10_elem_lapack_cpu const& _elem){

      if(__elem != NULL && __elemNum != 0)
        uni10_elem_free(__elem, __elemNum*sizeof(uni10_type));

      __elemNum =  _isdiag ? std::min(_Rnum, _Cnum) : _Rnum * _Cnum ;

      uni10_uint64 memsize = __elemNum * sizeof(uni10_type);

      if ( memsize ){

        __elem = (uni10_type*)uni10_elem_alloc( memsize );
        uni10_elem_copy( __elem, _elem.__elem, memsize );

      }

    };

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, const uni10_type* src){

      if(__elem != NULL && __elemNum != 0)
        uni10_elem_free(__elem, __elemNum*sizeof(uni10_type));

      __elemNum =  _isdiag ? std::min(_Rnum, _Cnum) : _Rnum * _Cnum ;

      uni10_uint64 memsize = __elemNum * sizeof(uni10_type);

      if ( memsize ){

        __elem = (uni10_type*)uni10_elem_alloc( memsize );
        if(src != NULL){
          uni10_elem_copy( __elem, src, memsize );
        }
        else{
          uni10_elemBzero( __elem, memsize );
        }

      }

    };

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::copy(uni10_elem_lapack_cpu const& _elem){

      uni10_uint64 memsize = __elemNum * sizeof(uni10_type);
      uni10_elem_copy( __elem, _elem.__elem, memsize );

    };

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::clear(){

      if(__elem != NULL && __elemNum != 0)
        uni10_elem_free(__elem, __elemNum * sizeof(uni10_type));

      __elemNum = 0;
      __elem = NULL;

    }

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::copy(uni10_uint64 begin_idx, const uni10_elem_lapack_cpu<uni10_type>& src, uni10_uint64 begin_src_idx, uni10_uint64 len){

      uni10_elem_copy(__elem + begin_idx, src.__elem + begin_src_idx, len*sizeof(uni10_type));

    }

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::copy(uni10_uint64 begin_idx, const uni10_elem_lapack_cpu<uni10_type>& src, uni10_uint64 len){

      this->copy(begin_idx, src, 0, len);

    }

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::catElem(const uni10_elem_lapack_cpu<uni10_type>& src){
        
      this->__elem = (uni10_type*)realloc(this->__elem, (this->__elemNum +src.__elemNum)*sizeof(uni10_type));
      this->copy(this->__elemNum, src, src.__elemNum);
      this->__elemNum += src.__elemNum;

    }

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::print_elem(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag) const{

      fprintf(stdout, "\n%ld x %ld = %ld [ Real ElemNum: %ld ]", _Rnum, _Cnum, _Rnum*_Cnum, __elemNum);

      if(__uni10_typeid == 1)  fprintf(stdout, ", REAL");
      else if(__uni10_typeid == 2)   fprintf(stdout, ", COMPLEX" );

      if(_isdiag)
        fprintf(stdout, ", Diagonal");

      fprintf(stdout, "\n\n");

      if ( _Rnum == 1 ) {
        fprintf(stdout, "[ " );
      }
      else {
        fprintf(stdout, "[\n" );
      }

      if(__elem == NULL){
        fprintf(stdout, "\nThe uni10_elem_lapack_cpu has not been allocated or linked. \n\n" );
        fprintf(stdout, "];\n" );
      }
      else if(_isdiag){
        for( uni10_int i = 0; i < (uni10_int)_Rnum; ++i ) {
          for( uni10_int j = 0; j < (uni10_int)_Cnum; ++j ) {
            if ( i != j) {
              if(__uni10_typeid == 2)
                fprintf(stdout, "   0.              " );
              else
                fprintf(stdout, "   0.    " );
            }
            else {
              uni10_print_elem_i(__elem[ i ]);
            }
          }
          if ( _Rnum > 1 ) 
            fprintf(stdout, "\n" );
          else 
            fprintf(stdout, " " );
        }
        fprintf(stdout, "];\n" );

      }
      else{
        for( uni10_int i = 0; i < (uni10_int)_Rnum; ++i ) {
          for( uni10_int j = 0; j < (uni10_int)_Cnum; ++j ) {
            if ( __elem[ i * _Cnum + j] == 0.) {
              if(__uni10_typeid == 2)
                fprintf(stdout, "   0.              " );
              else
                fprintf(stdout, "   0.    " );
            }
            else {
              uni10_print_elem_i(__elem[ i * _Cnum + j ]);
            }
          }
          if ( _Rnum > 1 ) 
            fprintf(stdout, "\n" );
          else 
            fprintf(stdout, " " );
        }
        fprintf(stdout, "];\n" );
      }

    }

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::save(FILE* fp) const{

      fwrite(&__uni10_typeid, sizeof(__uni10_typeid), 1, fp);
      fwrite(&__ongpu, sizeof(__ongpu), 1, fp);
      fwrite(&__elemNum, sizeof(__elemNum), 1, fp);
      fwrite(__elem, sizeof(uni10_type), __elemNum, fp);

    }

  //
  // This function has a requerment.
  // If the length of elements in the file is equal to the __elemNum, we copy the values directly without allocation because the operation malloc() will change the address.
  // If the address of __elem is changed, the UniTensor<T>::load() will crash.
  //
  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::load(FILE* fp){

      uni10_type_id buftype;
      uni10_error_msg(!fread(&buftype, sizeof(__uni10_typeid), 1, fp), "%s", "Loading __uni10_typeid is failure. (UNI10_LAPACKE_CPU<T>)");

      uni10_error_msg(buftype != __uni10_typeid, "%s", "TYPE ERROR. Can't loading a Real or Complex container to a Complex or Real one respectively.");

      uni10_error_msg(!fread(&__ongpu, sizeof(__ongpu), 1, fp), "%s", "Loading __ongpu is failure. (UNI10_LAPACKE_CPU<T>)");
      uni10_uint64 bufelemNum;
      uni10_error_msg(!fread(&bufelemNum, sizeof(__elemNum), 1, fp), "%s", "Loading __elemNum is failure. (UNI10_LAPACKE_CPU<T>)");

      if(__elem != NULL && __elemNum != 0 && bufelemNum != __elemNum)
        uni10_elem_free(__elem, __elemNum*sizeof(uni10_type));

      uni10_uint64 memsize = bufelemNum * sizeof(uni10_type);

      if ( memsize ){

        if(bufelemNum != __elemNum)
          __elem = (uni10_type*)uni10_elem_alloc( memsize );

        __elemNum = bufelemNum;
        uni10_error_msg(!fread(__elem, sizeof(uni10_type), __elemNum, fp), "%s", "Loading __elem is failure. (UNI10_LAPACKE_CPU<T>)");

      }

    }

  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::reshape(const std::vector<int>& ori_bdDims, const std::vector<int>& new_bdIdx, uni10_elem_lapack_cpu<uni10_type>& out_elem, bool inorder){

      int bondNum = ori_bdDims.size();

      if(!inorder){

        std::vector<int> newAcc(bondNum);
        newAcc[bondNum - 1] = 1;

        std::vector<int> transAcc(bondNum);
        transAcc[bondNum - 1] = 1;
        for(int b = bondNum - 1; b > 0; b--){
          newAcc[b - 1] = newAcc[b] * ori_bdDims[new_bdIdx[b]];
        }

        std::vector<int> newbondDims(bondNum);
        std::vector<int> idxs(bondNum);

        for(int b = 0; b < bondNum; b++){
          transAcc[new_bdIdx[b]] = newAcc[b];
          newbondDims[b] = ori_bdDims[new_bdIdx[b]];
          idxs[b] = 0;
        }

        int cnt_ot = 0;
        for(uni10_uint64 i = 0; i < __elemNum; i++){
          out_elem.__elem[cnt_ot] = this->__elem[i];
          for(int bend = bondNum - 1; bend >= 0; bend--){
            idxs[bend]++;
            if(idxs[bend] < ori_bdDims[bend]){
              cnt_ot += transAcc[bend];
              break;
            }
            else{
              cnt_ot -= transAcc[bend] * (idxs[bend] - 1);
              idxs[bend] = 0;
            }
          }
        }
      }else{
        out_elem.copy(0, *this, this->__elemNum);
      }

    }


  /*
  template<typename uni10_type>
    void uni10_elem_lapack_cpu<uni10_type>::resize(uni10_uint64 _row, uni10_uint64 _col, uni10_uint64& Rnum, uni10_uint64& Cnum, uni10_bool& isdiag, uni10_const_bool& _fixHead){

      if(_fixHead){

        if(isdiag){

          uni10_uint64 _elemNum = _row < _col ? _row : _col;

          if(_elemNum > __elemNum){

            uni10_type* _elem = (uni10_type*)uni10_elem_alloc(_elemNum * sizeof(uni10_type));

            uni10_elemBzero(_elem, _elemNum * sizeof(uni10_type));

            uni10_elem_copy(_elem, __elem, __elemNum * sizeof(uni10_type));

            //if(__elem != NULL && __elemNum == 0)
            if(__elem != NULL && __elemNum != 0)
              uni10_elem_free( __elem, __elemNum * sizeof(uni10_type) );

            __elem = _elem;

          }
          else
            shrinkWithoutFree( (__elemNum - _elemNum) * sizeof(uni10_type) );

          __elemNum = _elemNum;

          Rnum = _row;
          Cnum = _col;

        }
        else{

          uni10_uint64 _elemNum = _row * _col;

          if(_col == Cnum){

            if(_row > Rnum){

              uni10_type* _elem = (uni10_type*)uni10_elem_alloc( _elemNum * sizeof(uni10_type) );

              uni10_elemBzero( _elem, _elemNum * sizeof(uni10_type) );

              uni10_elem_copy(_elem, __elem, __elemNum * sizeof(uni10_type) );

              if(__elem != NULL)
                uni10_elem_free( __elem, __elemNum * sizeof(uni10_type) );

              __elem = _elem;

            }
            else
              shrinkWithoutFree( (__elemNum - _elemNum) * sizeof(uni10_type) );

            __elemNum = _elemNum;

            Rnum = _row;
          }
          else{

            uni10_uint64 data_row = _row < Rnum ? _row : Rnum;
            uni10_uint64 data_col = _col < Cnum ? _col : Cnum;

            uni10_type* _elem = (uni10_type*)uni10_elem_alloc( _elemNum * sizeof(uni10_type) );

            uni10_elemBzero( _elem, _elemNum * sizeof(uni10_type) );

            for(int r = 0; r < (int)data_row; r++)
              uni10_elem_copy( &(_elem[r * _col]), &(__elem[r * Cnum]), data_col * sizeof(uni10_type) );

            if(__elem != NULL)
              uni10_elem_free(__elem, __elemNum * sizeof(uni10_type) );

            __elem    = _elem;
            __elemNum = _elemNum;

            Rnum = _row;
            Cnum = _col;

          }

        }

      }else{

        uni10_error_msg(true, "%s", "Resize fixTail is developping !!!");

      }

    }
    */
  template<> template<>
    void uni10_elem_lapack_cpu<uni10_complex128>::init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, uni10_elem_lapack_cpu<uni10_double64> const& _elem){

      if(__elem != NULL && __elemNum != 0)
        uni10_elem_free(__elem, __elemNum*sizeof(uni10_complex128));

      __elemNum =  _isdiag ? std::min(_Rnum, _Cnum) : _Rnum * _Cnum ;

      uni10_uint64 memsize = __elemNum * sizeof(uni10_complex128);

      if ( memsize ){

        __elem = (uni10_complex128*)uni10_elem_alloc( memsize );
        uni10_elem_cast( __elem, _elem.__elem, __elemNum );

      }

    };

  template<> template<>
    void uni10_elem_lapack_cpu<uni10_double64>::init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, uni10_elem_lapack_cpu<uni10_complex128> const& _elem){

      if(__elem != NULL && __elemNum != 0)
        uni10_elem_free(__elem, __elemNum*sizeof(uni10_double64));

      __elemNum =  _isdiag ? std::min(_Rnum, _Cnum) : _Rnum * _Cnum ;

      uni10_uint64 memsize = __elemNum * sizeof(uni10_double64);

      if ( memsize ){

        __elem = (uni10_double64*)uni10_elem_alloc( memsize );
        uni10_elem_cast( __elem, _elem.__elem, __elemNum);

      }

    };


  template<> template<>
    uni10_elem_lapack_cpu<uni10_complex128>::uni10_elem_lapack_cpu(uni10_elem_lapack_cpu<uni10_double64> const& _elem):
      __uni10_typeid(_elem.__uni10_typeid), __ongpu(_elem.__ongpu), __elemNum(_elem.__elemNum), __elem(NULL){

      uni10_error_msg(_elem.__elem != NULL && _elem.__elemNum==0, "%s", "The pointer of this element container point to another pointer and does not allocate. Hence we can't copy it.");
      init(1, _elem.__elemNum, false, _elem);

    };

  template<> template<>
    uni10_elem_lapack_cpu<uni10_double64>::uni10_elem_lapack_cpu(uni10_elem_lapack_cpu<uni10_complex128> const& _elem):
      __uni10_typeid(_elem.__uni10_typeid), __ongpu(_elem.__ongpu), __elemNum(_elem.__elemNum), __elem(NULL){

      uni10_error_msg(_elem.__elem != NULL && _elem.__elemNum==0, "%s", "The pointer of this element container point to another pointer and does not allocate. Hence we can't copy it.");
      init(1, _elem.__elemNum, false,_elem);

    };

  template<> template<>
    uni10_elem_lapack_cpu<uni10_complex128>& uni10_elem_lapack_cpu<uni10_complex128>::operator=(uni10_elem_lapack_cpu<uni10_double64> const& _elem){

      uni10_error_msg(_elem.__elem != NULL && _elem.__elemNum==0, "%s", "The pointer of this element container point to another pointer and does not allocate. Hence we can't copy it.");
      __uni10_typeid = _elem.__uni10_typeid;
      __ongpu        = _elem.__ongpu;
      init(1, _elem.__elemNum, false, _elem);
      return *this;

    };

  template<> template<>
    uni10_elem_lapack_cpu<uni10_double64>& uni10_elem_lapack_cpu<uni10_double64>::operator=(uni10_elem_lapack_cpu<uni10_complex128> const& _elem){

      uni10_error_msg(_elem.__elem != NULL && _elem.__elemNum==0, "%s", "The pointer of this element container point to another pointer and does not allocate. Hence we can't copy it.");
      __uni10_typeid = _elem.__uni10_typeid;
      __ongpu        = _elem.__ongpu;
      init(1, _elem.__elemNum, false, _elem);
      return *this;

    };

  template<> template<>
    void uni10_elem_lapack_cpu<uni10_double64>::copy(uni10_elem_lapack_cpu<uni10_complex128> const& _elem){

      uni10_elem_cast( __elem, _elem.__elem, __elemNum);

    };

  template<> template<>
    void uni10_elem_lapack_cpu<uni10_complex128>::copy(uni10_elem_lapack_cpu<uni10_double64> const& _elem){

      uni10_elem_cast( __elem, _elem.__elem, __elemNum);

    };


  template class uni10_elem_lapack_cpu<uni10_double64>;
  template class uni10_elem_lapack_cpu<uni10_complex128>;

} /* namespace uni10 */
