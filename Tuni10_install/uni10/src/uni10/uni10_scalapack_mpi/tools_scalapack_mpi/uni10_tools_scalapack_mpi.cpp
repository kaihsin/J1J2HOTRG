#include "uni10/uni10_env_info.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_tools_scalapack_mpi.h"

namespace uni10{

  void* uni10_elem_alloc(uni10_uint64 memsize){

    void* ptr = NULL;
    ptr = malloc(memsize);

    uni10_error_msg(ptr==NULL, "%s","Fails in allocating memory.");
    return ptr;

  }

  void* uni10_elem_copy(void* des, const void* src, uni10_uint64 memsize){

    return memcpy(des, src, memsize);

  }

  void uni10_elem_free(void* ptr, uni10_uint64 memsize){

    free(ptr);

  }

  void uni10_elemBzero(void* ptr, uni10_uint64 memsize){

    memset(ptr, 0, memsize);

  }

  // For double 
  //
  void uni10_setDiag(uni10_double64* elem, uni10_double64* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }

  void uni10_getDiag(uni10_double64* elem, uni10_double64* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }

  void uni10_getUpTri(uni10_double64* elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }

  void uni10_getDnTri(uni10_double64* elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }

  void uni10_print_elem_i(const uni10_double64& elem_i){

    fprintf(stdout, " %8.4f", elem_i);

  }

  uni10_double64 UNI10_REAL( uni10_double64 elem_i ){

    return elem_i;

  }

  void uni10_getUpTri(uni10_complex128* elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }

  void uni10_getDnTri(uni10_complex128* elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }

  uni10_double64 UNI10_IMAG( uni10_double64 elem_i ){

    return 0.;

  }

  // For complex 
  //
  void uni10_setDiag(uni10_complex128* elem, uni10_complex128* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }
  void uni10_getDiag(uni10_complex128* elem, uni10_complex128* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }

  void uni10_print_elem_i(const uni10_complex128& elem_i){

    fprintf(stdout, " %8.4f+%8.4fi", Z_REAL( elem_i ), Z_IMAG( elem_i ) );

  }

  uni10_double64 UNI10_REAL( uni10_complex128 elem_i ){

    return elem_i.real();

  }

  uni10_double64 UNI10_IMAG( uni10_complex128 elem_i ){

    return elem_i.imag();

  }

  void ToReal(uni10_double64& M_i, uni10_double64 val){

    M_i = val;

  }

  void ToReal(uni10_complex128& M_i, uni10_double64 val){

    M_i.real(val);

  }

  void ToComplex(uni10_double64& M_i, uni10_double64 val){

    // Do nothing

  }

  void ToComplex(uni10_complex128& M_i, uni10_double64 val){

    M_i.imag(val);

  }

  // Convert
  void uni10_elem_cast(uni10_complex128* des, uni10_double64* src, uni10_uint64 N){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }

  void uni10_elem_cast(uni10_double64* des, uni10_complex128* src, uni10_uint64 N){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }

  void shrinkWithoutFree(uni10_uint64 memsize){

    uni10_error_msg(true, "%s", "Developping!!!\n");

  }

  void printDesc(const uni10_int* desc){

    for(uni10_int i = 0; i < env_variables.get_info().nprocs_mpi; i++ ){

      if(env_variables.get_info().rank_mpi == i){
        printf("\n/********** Rank %d **********\n", env_variables.get_info().rank_mpi );
        printf("/* (GLOBAL) DTYPE : %d\n", desc[0]);
        printf("/* (GLOBAL) CTXT  : %d\n", desc[1]);
        printf("/* (GLOBAL) M     : %d\n", desc[2]);
        printf("/* (GLOBAL) N     : %d\n", desc[3]);
        printf("/* (GLOBAL) MB    : %d\n", desc[4]);
        printf("/* (GLOBAL) NB    : %d\n", desc[5]);
        printf("/* (GLOBAL) RSRC  : %d\n", desc[6]);
        printf("/* (GLOBAL) CSRC  : %d\n", desc[7]);
        printf("/* (LODCAL) LLD   : %d\n", desc[8]);
        printf("/*****************************\n");
      }

      MPI_Barrier(MPI_COMM_WORLD);

    }

  }

  void broadcast(uni10_double64* dest, const uni10_double64* buffer, const uni10_int count, const uni10_int root){

    if(env_variables.get_info().rank_mpi == root)
      memcpy(dest, buffer, count*sizeof(uni10_double64));

    MPI_Bcast(dest, count, MPI_DOUBLE, root, MPI_COMM_WORLD);

  };

  void broadcast(uni10_complex128* dest, const uni10_complex128* buffer, const uni10_int count, const uni10_int root){

    if(env_variables.get_info().rank_mpi == root)
      memcpy(dest, buffer, count*sizeof(uni10_complex128));

    MPI_Bcast(dest, count, MPI_C_DOUBLE_COMPLEX, root, MPI_COMM_WORLD);

  };

  void split_mast2dist(const uni10_double64* src, const uni10_int* r_lens, const uni10_int& nprow, const uni10_int* c_lens, const uni10_int& npcol, 
      const uni10_int* gdists, const uni10_int& Rnum_, const uni10_int& src_root, 
      uni10_double64* dist, const uni10_int& blkrow, const uni10_int& blkcol, const uni10_int& rank){

    MPI_Request request_s;
    MPI_Request request_r;
    MPI_Status status;

    if(rank==src_root){

      uni10_error_msg( src  == NULL, "%s", "The source ptr is NULL.");

      uni10_int k = 0;

      for(uni10_int i = 0; i < npcol; i++){

        for(uni10_int j = 0; j < nprow; j++){

          MPI_Datatype blocktype;
          MPI_Datatype blocktype2;

          MPI_Type_vector(c_lens[i], r_lens[j], Rnum_, MPI_DOUBLE, &blocktype2);
          MPI_Type_create_resized( blocktype2, 0, sizeof(MPI_DOUBLE), &blocktype );
          MPI_Type_commit(&blocktype);
          MPI_Isend(src+gdists[k], 1, blocktype, k, k, MPI_COMM_WORLD, &request_s);
          k++;

        }

      }

    }

    MPI_Irecv(dist, blkrow*blkcol, MPI_DOUBLE, src_root, rank, MPI_COMM_WORLD, &request_r);
    MPI_Wait(&request_r, &status);

  }


  void split_mast2dist(const uni10_complex128* src, const uni10_int* r_lens, const uni10_int& nprow, const uni10_int* c_lens, const uni10_int& npcol, 
      const uni10_int* gdists, const uni10_int& Rnum_, const uni10_int& src_root, 
      uni10_complex128* dist, const uni10_int& blkrow, const uni10_int& blkcol, const uni10_int& rank){

    MPI_Request request_s;
    MPI_Request request_r;
    MPI_Status status;

    if(rank==src_root){

      uni10_error_msg( src  == NULL, "%s", "The source ptr is NULL.");

      uni10_int k = 0;

      for(uni10_int i = 0; i < npcol; i++)

        for(uni10_int j = 0; j < nprow; j++){

          MPI_Datatype blocktype;
          MPI_Datatype blocktype2;

          MPI_Type_vector(c_lens[i], r_lens[j], Rnum_, MPI_C_DOUBLE_COMPLEX, &blocktype2);
          MPI_Type_create_resized( blocktype2, 0, sizeof(MPI_C_DOUBLE_COMPLEX), &blocktype );
          MPI_Type_commit(&blocktype);
          MPI_Isend(src+gdists[k], 1, blocktype, k, k, MPI_COMM_WORLD, &request_s);
          k++;

        }

    }

    MPI_Irecv(dist, blkrow*blkcol, MPI_C_DOUBLE_COMPLEX, src_root, rank, MPI_COMM_WORLD, &request_r);
    MPI_Wait(&request_r, &status);

  }

  void gather_dist2mast(const uni10_double64* dist, const uni10_int& blkrow, const uni10_int& blkcol, const uni10_int& rank, 
      uni10_double64* dest_ptr, const uni10_int* r_lens, const uni10_int& nprow, const uni10_int* const c_lens, const uni10_int& npcol, 
      const uni10_int* gdists, const uni10_int& Rnum_, const uni10_int& dest){

    MPI_Request request_s;
    MPI_Request request_r;
    MPI_Status status;

    MPI_Isend(dist, blkrow*blkcol, MPI_DOUBLE, dest, rank, MPI_COMM_WORLD, &request_s);

    if(rank==dest){

      uni10_int k = 0;

      for(uni10_int i = 0; i < npcol; i++){

        for(uni10_int j = 0; j < nprow; j++){

          MPI_Datatype blocktype;
          MPI_Datatype blocktype2;
          MPI_Type_vector(c_lens[i], r_lens[j], Rnum_, MPI_DOUBLE, &blocktype2);
          MPI_Type_create_resized( blocktype2, 0, sizeof(MPI_DOUBLE), &blocktype );
          MPI_Type_commit(&blocktype);
          MPI_Irecv(dest_ptr+gdists[k], 1, blocktype, k, k, MPI_COMM_WORLD, &request_r);
          MPI_Wait(&request_r, &status);
          k++;
        }

      }

    }

  }

  void gather_dist2mast(const uni10_complex128* dist, const uni10_int& blkrow, const uni10_int& blkcol, const uni10_int& rank, 
      uni10_complex128* dest_ptr, const uni10_int* r_lens, const uni10_int& nprow, const uni10_int* c_lens, const uni10_int& npcol, 
      const uni10_int* gdists, const uni10_int& Rnum_, const uni10_int& dest){

    MPI_Request request_s;
    MPI_Request request_r;
    MPI_Status status;

    MPI_Isend(dist, blkrow*blkcol, MPI_C_DOUBLE_COMPLEX, dest, rank, MPI_COMM_WORLD, &request_s);

    if(rank==dest){

      uni10_int k = 0;

      for(uni10_int i = 0; i < npcol; i++){

        for(uni10_int j = 0; j < nprow; j++){

          MPI_Datatype blocktype;
          MPI_Datatype blocktype2;
          MPI_Type_vector(c_lens[i], r_lens[j], Rnum_, MPI_C_DOUBLE_COMPLEX, &blocktype2);
          MPI_Type_create_resized( blocktype2, 0, sizeof(MPI_C_DOUBLE_COMPLEX), &blocktype );
          MPI_Type_commit(&blocktype);
          MPI_Irecv(dest_ptr+gdists[k], 1, blocktype, k, k, MPI_COMM_WORLD, &request_r);
          MPI_Wait(&request_r, &status);
          k++;
        }

      }

    }

  }

} /* namespace uni10 */
