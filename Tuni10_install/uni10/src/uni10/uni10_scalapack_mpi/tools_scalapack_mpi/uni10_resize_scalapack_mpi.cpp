#include "uni10/uni10_type.h"
#include "uni10/uni10_scalapack_mpi/tools_scalapack_mpi/uni10_resize_scalapack_mpi.h"

namespace uni10{

  void resize_(uni10_elem_double64& Eout, const uni10_elem_double64& Ein, uni10_const_bool& _fixHead){

    if(_fixHead){

      uni10_double64* global_elem = NULL;
      if(Ein.rank==Ein.master)
        global_elem = (uni10_double64*)malloc(Eout.Rnum_*Eout.Cnum_*sizeof(uni10_double64));

      ////Gather
      MPI_Request request_s;
      MPI_Request request_r;
      MPI_Status status_mpi;

      uni10_int blkrow = Ein.r_offset <= (uni10_int)Eout.Rnum_ ? Ein.blockrow : Eout.Rnum_ - Ein.r_head;
      uni10_int blkcol = Ein.c_offset <= (uni10_int)Eout.Cnum_ ? Ein.blockcol : Eout.Cnum_ - Ein.c_head;
      blkrow = blkrow < 0 ? 0 : blkrow;
      blkcol = blkcol < 0 ? 0 : blkcol;

      for(uni10_int i = 0; i < Ein.npnum; i++){
        if(Ein.rank==i){
          //    //printf("S: rank: %d, blkrow: %d, blkcol: %d\n", rank, blkrow, blkcol);
          MPI_Datatype blocktype;
          MPI_Datatype blocktype2;
          MPI_Type_vector(blkcol, blkrow, Ein.blockrow, MPI_DOUBLE, &blocktype2);
          MPI_Type_create_resized( blocktype2, 0, sizeof(MPI_DOUBLE), &blocktype );
          MPI_Type_commit(&blocktype);
          MPI_Isend(Ein.__elem, 1, blocktype, Ein.master, Ein.rank, MPI_COMM_WORLD, &request_s);
        }
      }

      if(Ein.rank==Ein.master){
        uni10_int c_end=0, c_end_old = 0, gdisp = 0, k=0;
        for(uni10_int c = 0; c < Ein.npcol; c++){
          c_end += Ein.c_lens[c];
          uni10_int r_end=0, r_end_old = 0;
          for(uni10_int r = 0; r < Ein.nprow; r++){
            MPI_Datatype blocktype;
            MPI_Datatype blocktype2;
            r_end += Ein.r_lens[r];
            blkrow = r_end <= (uni10_int)Eout.Rnum_ ? Ein.r_lens[r] : Eout.Rnum_ - r_end_old;
            blkcol = c_end <= (uni10_int)Eout.Cnum_ ? Ein.c_lens[c] : Eout.Cnum_ - c_end_old;
            blkrow = blkrow < 0 ? 0 : blkrow;
            blkcol = blkcol < 0 ? 0 : blkcol;
            gdisp  = c_end_old*Eout.Rnum_+r_end_old;
            MPI_Type_vector(blkcol, blkrow, Eout.Rnum_, MPI_DOUBLE, &blocktype2);
            MPI_Type_create_resized( blocktype2, 0, sizeof(MPI_DOUBLE), &blocktype );
            MPI_Type_commit(&blocktype);
            MPI_Irecv(global_elem+gdisp, 1, blocktype , k, k, MPI_COMM_WORLD, &request_r);
            MPI_Wait(&request_r, &status_mpi);
            r_end_old = r_end;
            k++;
          }
          c_end_old = c_end;
        }

      }

      Eout.setElem(global_elem);
      free(global_elem);

    }else{

      uni10_error_msg(true, "%s", "Resize fixTail is developping !!!");

    }

  }

  void resize_(uni10_elem_complex128& Eout, const uni10_elem_complex128& Ein, uni10_const_bool& _fixHead){

    if(_fixHead){

      uni10_complex128* global_elem = NULL;
      if(Ein.rank==Ein.master)
        global_elem = (uni10_complex128*)malloc(Eout.Rnum_*Eout.Cnum_*sizeof(uni10_complex128));

      ////Gather
      MPI_Request request_s;
      MPI_Request request_r;
      MPI_Status status_mpi;

      uni10_int blkrow = Ein.r_offset <= (uni10_int)Eout.Rnum_ ? Ein.blockrow : Eout.Rnum_ - Ein.r_head;
      uni10_int blkcol = Ein.c_offset <= (uni10_int)Eout.Cnum_ ? Ein.blockcol : Eout.Cnum_ - Ein.c_head;
      blkrow = blkrow < 0 ? 0 : blkrow;
      blkcol = blkcol < 0 ? 0 : blkcol;

      for(uni10_int i = 0; i < Ein.npnum; i++){
        if(Ein.rank==i){
          //    //printf("S: rank: %d, blkrow: %d, blkcol: %d\n", rank, blkrow, blkcol);
          MPI_Datatype blocktype;
          MPI_Datatype blocktype2;
          MPI_Type_vector(blkcol, blkrow, Ein.blockrow, MPI_C_DOUBLE_COMPLEX, &blocktype2);
          MPI_Type_create_resized( blocktype2, 0, sizeof(MPI_C_DOUBLE_COMPLEX), &blocktype );
          MPI_Type_commit(&blocktype);
          MPI_Isend(Ein.__elem, 1, blocktype, Ein.master, Ein.rank, MPI_COMM_WORLD, &request_s);
        }
      }

      if(Ein.rank==Ein.master){
        uni10_int c_end=0, c_end_old = 0, gdisp = 0, k=0;
        for(uni10_int c = 0; c < Ein.npcol; c++){
          c_end += Ein.c_lens[c];
          uni10_int r_end=0, r_end_old = 0;
          for(uni10_int r = 0; r < Ein.nprow; r++){
            MPI_Datatype blocktype;
            MPI_Datatype blocktype2;
            r_end += Ein.r_lens[r];
            blkrow = r_end <= (uni10_int)Eout.Rnum_ ? Ein.r_lens[r] : Eout.Rnum_ - r_end_old;
            blkcol = c_end <= (uni10_int)Eout.Cnum_ ? Ein.c_lens[c] : Eout.Cnum_ - c_end_old;
            blkrow = blkrow < 0 ? 0 : blkrow;
            blkcol = blkcol < 0 ? 0 : blkcol;
            gdisp  = c_end_old*Eout.Rnum_+r_end_old;
            MPI_Type_vector(blkcol, blkrow, Eout.Rnum_, MPI_C_DOUBLE_COMPLEX, &blocktype2);
            MPI_Type_create_resized( blocktype2, 0, sizeof(MPI_C_DOUBLE_COMPLEX), &blocktype );
            MPI_Type_commit(&blocktype);
            MPI_Irecv(global_elem+gdisp, 1, blocktype , k, k, MPI_COMM_WORLD, &request_r);
            MPI_Wait(&request_r, &status_mpi);
            r_end_old = r_end;
            k++;
          }
          c_end_old = c_end;
        }

      }

      Eout.setElem(global_elem);
      free(global_elem);

    }else{

      uni10_error_msg(true, "%s", "Resize fixTail is developping !!!");

    }

  }

  ////  for(int i = 0; i < _row; i++){
  ////    for(int j = 0; j < _col; j++)
  ////      std::cout <<  global_elem[i*_col+j] << " ";
  ////    std::cout << std::endl;
  ////  }

  ////  MPI_Barrier(MPI_COMM_WORLD);


}
