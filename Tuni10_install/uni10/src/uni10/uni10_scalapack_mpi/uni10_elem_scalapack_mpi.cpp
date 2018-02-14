#include "uni10/uni10_type.h"
#include "uni10/uni10_error.h"
#include "uni10/uni10_scalapack_mpi/uni10_elem_scalapack_mpi.h"

namespace uni10{

  // Done
  template<typename uni10_type>
    uni10_elem_scalapack_mpi<uni10_type>::uni10_elem_scalapack_mpi(): __uni10_typeid(UNI10_TYPE_ID(uni10_type)), __ongpu(false), __elemNum(0), __elem(NULL), status(0){
      
    };

  // Done
  template<typename uni10_type>
    uni10_elem_scalapack_mpi<uni10_type>::uni10_elem_scalapack_mpi(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, uni10_bool _ongpu): __uni10_typeid(UNI10_TYPE_ID(uni10_type)), __ongpu(_ongpu),__elem(NULL), status(0){

      init(_Rnum, _Cnum, _isdiag, NULL);

    }

  // Done
  template<typename uni10_type>
    uni10_elem_scalapack_mpi<uni10_type>::uni10_elem_scalapack_mpi(const uni10_type* src, uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, uni10_bool _ongpu): __uni10_typeid(UNI10_TYPE_ID(uni10_type)), __ongpu(_ongpu), __elem(NULL), status(0){

      init(_Rnum, _Cnum, _isdiag, src);

    };

  // Done
  template<typename uni10_type>
    uni10_elem_scalapack_mpi<uni10_type>::uni10_elem_scalapack_mpi(const uni10_elem_scalapack_mpi& _elem){

      status = 0;
      this->copy(_elem);

    };

  // Done
  template<typename uni10_type>
    uni10_elem_scalapack_mpi<uni10_type>::~uni10_elem_scalapack_mpi(){

      bool _isdiag = (uni10_int)__elemNum != Rnum_ * Cnum_;

      if(__elem != NULL && __elemNum != 0)
        uni10_elem_free(__elem, __elemNum * sizeof(uni10_type));

      __elem    = NULL;
      __elemNum = 0;

      if(status && rank==master && !_isdiag){
        free(r_lens); free(c_lens);
        free(rdists); free(cdists); free(gdists);
      }

      MPI_Barrier(MPI_COMM_WORLD);

    };

  // Done
  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::set_zeros(){

      uni10_error_msg( status==0, "%s", "Please initialize the uni10_elem with the constructor uni10(uni10_uint64, uni10_uint64, bool) befero setting the elements.");

      uni10_uint64 memsize = blockrow*blockcol* sizeof(uni10_type);

      uni10_elemBzero( __elem, memsize );

    };

  // Done
  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::setElem(const uni10_type* src, bool src_dist){

      uni10_error_msg( status == 0, "%s", "Please initialize the uni10_elem with the constructor uni10(uni10_uint64, uni10_uint64, bool) befero setting the elements.");
      if(rank==master)
        uni10_error_msg( src  == NULL, "%s", "The source ptr is NULL.");

      bool _isdiag = (uni10_int)__elemNum != Rnum_ * Cnum_;

      if(!_isdiag){

        if(!src_dist)
          split_mast2dist(src, r_lens, nprow, c_lens, npcol, gdists, Rnum_, master, __elem, blockrow,blockcol, rank);
        else{
          uni10_elem_copy( __elem, src, blockrow*blockcol * sizeof(uni10_type) );
        }

      }else{

        if(!src_dist)
          broadcast(__elem, src, __elemNum, master);
        else{
          uni10_elem_copy( __elem, src, __elemNum * sizeof(uni10_type) );
        }

      }

    };

  // Done
  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::init(uni10_uint64 _Rnum, uni10_uint64 _Cnum, uni10_bool _isdiag, const uni10_type* src){

      if(status)
        this->clear();

      status = (_Rnum*_Cnum != 0) ? 1 : 0; // Initialization.

      // Total element number.
      __elemNum =  _isdiag ? std::min(_Rnum, _Cnum) : _Rnum * _Cnum ;

      master = env_variables.get_info().master;

      // C2f
      Rnum_ =  _Cnum;
      Cnum_ =  _Rnum;

      rank     = env_variables.get_info().rank_mpi;
      npnum    = env_variables.get_info().nprocs_mpi;

      ictxt    = env_variables.get_info().ictxt;

      nprow    = env_variables.get_info().nprow;
      npcol    = env_variables.get_info().npcol;

      myrow     = env_variables.get_info().myrow;
      mycol     = env_variables.get_info().mycol;

      rgrid    = ( env_variables.get_info().blockgrid > Rnum_) ? Rnum_ : env_variables.get_info().blockgrid; 
      cgrid    = ( env_variables.get_info().blockgrid > Cnum_) ? Cnum_ : env_variables.get_info().blockgrid; 

      uni10_error_msg(!(myrow>-1)&(mycol>-1)&(myrow<nprow)&(mycol<npcol), "%s", "Fail to initialize blocks' grids.");

      uni10_int memsize = 0;

      // If Matrix is diagonal, we put elements on global memory.
      // Need not  descriptor.
      if(!_isdiag){

        uni10_int izero = 0;

        blockrow    = numroc( &Rnum_   , &rgrid, &myrow, &izero, &nprow );
        blockcol    = numroc( &Cnum_   , &cgrid, &mycol, &izero, &npcol );

        memsize = blockrow*blockcol * sizeof(uni10_type);

        uni10_int itemp = std::max(1, blockrow);
        uni10_int info;

        // Generate a matrix descriptor.
        descinit( desc, &Rnum_ , &Cnum_  , &rgrid,  &cgrid, &izero, &izero, &ictxt, &itemp, &info );

        this->init_dists();

      }else{

        blockrow = 1;
        blockcol = __elemNum;
        memsize  = __elemNum * sizeof(uni10_type);

      }

      if ( memsize ){

        __elem = (uni10_type*)uni10_elem_alloc( memsize );

        if(src != NULL){
          uni10_elem_copy( __elem, src, memsize );
        }
        else{
          uni10_elemBzero( __elem, memsize );
        }

      }

      //printf("rank: %d, rgrid: %d, cgrid: %d, myrow: %d, mycol: %d, nprow: %d, npcol: %d, blockrow: %d, blockcol: %d\n", 
      //    rank, rgrid, cgrid, myrow, mycol, nprow, npcol, blockrow, blockcol);
      MPI_Barrier(MPI_COMM_WORLD);
      //exit(0);

    };

  // Done
  template <typename uni10_type> 
    void uni10_elem_scalapack_mpi<uni10_type>::copy(const uni10_elem_scalapack_mpi& _elem){

      if(status)
        this->clear();

      bool _isdiag = _elem.Rnum_ * _elem.Cnum_ != _elem.__elemNum;

      __uni10_typeid  = _elem.__uni10_typeid;
      __ongpu         = false;
      __elemNum       = _elem.__elemNum;

      ictxt     = _elem.ictxt;
      rank      = _elem.rank;
      npnum     = _elem.npnum;
      nprow     = _elem.nprow;
      npcol     = _elem.npcol;
      blockrow  = _elem.blockrow;
      blockcol  = _elem.blockcol;
      rgrid     = _elem.rgrid; 
      cgrid     = _elem.cgrid;
      myrow     = _elem.myrow;
      mycol     = _elem.mycol;
      Rnum_     = _elem.Rnum_;
      Cnum_     = _elem.Cnum_;
      master    = _elem.master;
      status    = _elem.status;

      uni10_error_msg(!(myrow>-1)&(mycol>-1)&(myrow<nprow)&(mycol<npcol), "%s", "Fail to initialize blocks' grids.");

      uni10_uint64 memsize = 0;

      if(!_isdiag){

        r_head    = _elem.r_head;
        r_offset  = _elem.r_offset;
        c_head    = _elem.c_head;
        c_offset  = _elem.c_offset;

        uni10_elem_copy(desc, _elem.desc, 9*sizeof(uni10_int));

        if(rank==master){

          memsize = nprow*sizeof(uni10_int);
          r_lens = (uni10_int*)malloc( memsize );
          uni10_elem_copy(r_lens, _elem.r_lens, memsize);

          memsize = npcol*sizeof(uni10_int);
          c_lens = (uni10_int*)malloc( memsize );
          uni10_elem_copy(c_lens, _elem.c_lens, memsize);

          memsize = nprow*sizeof(uni10_int);
          rdists = (uni10_int*)malloc( memsize );
          uni10_elem_copy(rdists, _elem.rdists, memsize);

          memsize = npcol*sizeof(uni10_int);
          cdists = (uni10_int*)malloc( memsize );
          uni10_elem_copy(cdists, _elem.cdists, memsize);

          memsize = npnum*sizeof(uni10_int);
          gdists = (uni10_int*)malloc( memsize );
          uni10_elem_copy(gdists, _elem.gdists, memsize);

        }

      }

      memsize = blockrow*blockcol * sizeof(uni10_type);

      if ( memsize ){

        __elem = (uni10_type*)uni10_elem_alloc( memsize );
        uni10_elem_copy( __elem, _elem.__elem, memsize );

      }

      MPI_Barrier(MPI_COMM_WORLD);

    };


  // Done
  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::init_dists(){

      if(rank==master){
        r_lens = (uni10_int*)malloc(nprow*sizeof(uni10_int));
        c_lens = (uni10_int*)malloc(npcol*sizeof(uni10_int));
      }

      MPI_Request request_sr;
      MPI_Request request_sc;
      MPI_Request request_rr;
      MPI_Request request_rc;
      MPI_Status  status_mpi, status_mpio ;

      if(rank<nprow)
        MPI_Isend(&blockrow, 1, MPI_INT, master, rank, MPI_COMM_WORLD, &request_sr);
      if(rank%nprow==0)
        MPI_Isend(&blockcol, 1, MPI_INT, master, rank, MPI_COMM_WORLD, &request_sc);

      if(rank==master){
        for(uni10_int r = 0; r < nprow; r++){
          MPI_Irecv(r_lens+r, 1, MPI_INT, r, r, MPI_COMM_WORLD, &request_rr);
          MPI_Wait(&request_rr, &status_mpi);
        }
        for(uni10_int c = 0; c < npcol; c++){
          uni10_int src_rank =  c*nprow;
          MPI_Irecv(c_lens+c, 1, MPI_INT, src_rank, src_rank, MPI_COMM_WORLD, &request_rc);
          MPI_Wait(&request_rc, &status_mpi);
        }
      }

      MPI_Request request_sro;
      MPI_Request request_rro;

      // Send r_head and r_offset
      if(rank==master){
        uni10_int p = 0;
        uni10_int head = 0; 
        uni10_int offset = 0; 
        for(uni10_int r = 0; r < nprow; r++){
          offset += r_lens[r];
          for(uni10_int c = 0; c < npcol; c++){
            p = c * nprow + r;
            if(p == master){
              r_head = head;
              r_offset = offset;
              continue;
            }
            MPI_Isend(&head  , 1, MPI_INT, p, 2*p  , MPI_COMM_WORLD, &request_sr);
            MPI_Isend(&offset, 1, MPI_INT, p, 2*p+1, MPI_COMM_WORLD, &request_sro);
          }
          head += r_lens[r];
        }
      }

      if(rank!=master){
        MPI_Irecv(&r_head  , 1, MPI_INT, master, 2*rank  , MPI_COMM_WORLD, &request_rr);
        MPI_Irecv(&r_offset, 1, MPI_INT, master, 2*rank+1, MPI_COMM_WORLD, &request_rro);
        MPI_Wait(&request_rr, &status_mpi);
        MPI_Wait(&request_rro, &status_mpio);
      }

      // Send c_head and c_offset
      if(rank==master){
        uni10_int p = 0;
        uni10_int head = 0; 
        uni10_int offset = 0; 
        for(uni10_int c = 0; c < npcol; c++){
          offset += c_lens[c];
          for(uni10_int r = 0; r < nprow; r++){
            p = c * nprow + r;
            if(p == master){
              c_head = head;
              c_offset = offset;
              continue;
            }
            MPI_Isend(&head  , 1, MPI_INT, p, 2*p  , MPI_COMM_WORLD, &request_sr);
            MPI_Isend(&offset, 1, MPI_INT, p, 2*p+1, MPI_COMM_WORLD, &request_sro);
          }
          head += c_lens[c];
        }
      }

      if(rank!=master){
        MPI_Irecv(&c_head  , 1, MPI_INT, master, 2*rank  , MPI_COMM_WORLD, &request_rr);
        MPI_Irecv(&c_offset, 1, MPI_INT, master, 2*rank+1, MPI_COMM_WORLD, &request_rro);
        MPI_Wait(&request_rr, &status_mpi);
        MPI_Wait(&request_rro, &status_mpio);
      }

      //printf("rank: %d, c_head: %d, c_offset: %d\n", rank, c_head, c_offset);
      //MPI_Barrier(MPI_COMM_WORLD);

      // Compute rdists, cdists and global dists
      if(rank==master){

        rdists = (uni10_int*)malloc(nprow*sizeof(uni10_int));
        cdists = (uni10_int*)malloc(npcol*sizeof(uni10_int));
        gdists = (uni10_int*)malloc(npnum*sizeof(uni10_int));

        // count rdists ( by row )
        rdists[0]=0;
        for(uni10_int r = 1; r < npcol; r++){
          rdists[r] = rdists[r-1] + c_lens[r-1];
          //printf("rdists[%d]: %d\n", r, rdists[r]);
        }

        // count cdists ( by row )
        cdists[0]=0;
        for(uni10_int c = 1; c < nprow; c++){
          cdists[c] = cdists[c-1] + r_lens[c-1];
          //printf("cdists[%d]: %d\n", c, cdists[c]);
        }

        uni10_int k = 0;
        for (uni10_int i=0; i<npcol; i++) {
          for (uni10_int j=0; j<nprow; j++) {
            gdists[k] = rdists[i]* Rnum_ + cdists[j];
            //printf("gdists[%d]: %d\n", k, gdists[k]);
            k++;
          }
        }

      }

    }

  // Done
  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::assign(uni10_uint64& _Rnum, uni10_uint64& _Cnum){

      this->init(_Rnum, _Cnum, false);

    }

  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::clear(){

      bool _isdiag = __elemNum != Rnum_ * Cnum_;

      if(__elem != NULL && __elemNum != 0)
        uni10_elem_free(__elem, __elemNum * sizeof(uni10_type));

      __elemNum = 0;
      __elem = NULL;


      if(status && rank==master && !_isdiag){
        free(r_lens); free(c_lens);
        free(rdists); free(cdists); free(gdists);
      }

      status = 0;

      MPI_Barrier(MPI_COMM_WORLD);

    }

  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::copy(uni10_uint64 begin_idx, const uni10_elem_scalapack_mpi<uni10_type>& src, uni10_uint64 begin_src_idx, uni10_uint64 len){

      uni10_elem_copy(__elem + begin_idx, src.__elem + begin_src_idx, len*sizeof(uni10_type));

    }

  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::copy(uni10_uint64 begin_idx, const uni10_elem_scalapack_mpi<uni10_type>& src, uni10_uint64 len){

      this->copy(begin_idx, src, 0, len);

    }

  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::catElem(const uni10_elem_scalapack_mpi<uni10_type>& src){

      this->__elem = (uni10_type*)realloc(this->__elem, (this->__elemNum +src.__elemNum)*sizeof(uni10_type));
      this->copy(this->__elemNum, src, src.__elemNum);
      this->__elemNum += src.__elemNum;

    }

  // Done
  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::print_elem(uni10_uint64 _Rnum, uni10_uint64 _Cnum, bool _isdiag) const{

      uni10_type* gather_elem = NULL;

      if(rank == master){
        gather_elem = (uni10_type*)uni10_elem_alloc(__elemNum*sizeof(uni10_type));
      }

      if(!_isdiag)
        gather_dist2mast(__elem, blockrow, blockcol, rank, gather_elem, r_lens, nprow, c_lens, npcol, gdists, Rnum_, master);
      else
        if(rank == master)
          uni10_elem_copy(gather_elem, __elem, __elemNum*sizeof(uni10_type));

      if(rank==master){

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

        if(gather_elem == NULL){
          fprintf(stdout, "\nThe uni10_elem_scalapack_mpi has not been allocated or linked. \n\n" );
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
                uni10_print_elem_i(gather_elem[ i ]);
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
              if ( gather_elem[ i * _Cnum + j] == 0.) {
                if(__uni10_typeid == 2)
                  fprintf(stdout, "   0.              " );
                else
                  fprintf(stdout, "   0.    " );
              }
              else {
                uni10_print_elem_i(gather_elem[ i * _Cnum + j ]);
              }
            }
            if ( _Rnum > 1 ) 
              fprintf(stdout, "\n" );
            else 
              fprintf(stdout, " " );
          }

          fprintf(stdout, "];\n" );
        }

        free(gather_elem);

      }

      MPI_Barrier(MPI_COMM_WORLD);

    }

  template<typename uni10_type>
    void uni10_elem_scalapack_mpi<uni10_type>::save(FILE* fp) const{

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
    void uni10_elem_scalapack_mpi<uni10_type>::load(FILE* fp){

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

  template class uni10_elem_scalapack_mpi<uni10_double64>;
  template class uni10_elem_scalapack_mpi<uni10_complex128>;

} /* namespace uni10 */
