#include "uni10/uni10_env_info/uni10_env_info_scalapack_mpi.h"

namespace uni10{


  uni10_env_mpi env_variables;

  void uni10_env_mpi::init(int argc, char** argv){

    // MPI paramters
    uni10_int& _rank_mpi   = uni10_sys_info.rank_mpi;
    uni10_int& _nprocs_mpi = uni10_sys_info.nprocs_mpi;
    uni10_int& _master     = uni10_sys_info.master;
    // BLACS parameter.
    uni10_int& _ictxt      = uni10_sys_info.ictxt;
    uni10_int& _nprow      = uni10_sys_info.nprow;
    uni10_int& _npcol      = uni10_sys_info.npcol;
    uni10_int& _blockgrid  = uni10_sys_info.blockgrid;
    uni10_int& _myrow      = uni10_sys_info.myrow;
    uni10_int& _mycol      = uni10_sys_info.mycol;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &_nprocs_mpi);

    if(!this->load_uni10_rc(_rank_mpi, _nprow, _npcol, _blockgrid, _master)){

      uni10_sys_info = sysinfo_mpi();

    }

    // Might be changed 
    uni10_int ConTxt = -1;
    uni10_int what = 0;
    //

    Cblacs_pinfo( &_rank_mpi, &_nprocs_mpi ) ;
    Cblacs_get( ConTxt, what, &_ictxt );
    Cblacs_gridinit( &_ictxt, (char*)"Col", _nprow, _npcol );
    // Get each blocks' myrow and mycol.
    Cblacs_gridinfo( _ictxt, &_nprow, &_npcol, &_myrow, &_mycol );
    
  }

  void uni10_env_mpi::clear(){

    etype         = mpi;

  }

  uni10_env_mpi::~uni10_env_mpi(){

    //uni10_sys_info.clear();
    Cblacs_barrier( uni10_sys_info.ictxt, (char*)"All" );
    Cblacs_gridexit( uni10_sys_info.ictxt );

    MPI_Finalize();

  }

  void uni10_env_mpi::used_memsize(const uni10_uint64& memsize){


  }

  bool uni10_env_mpi::load_uni10_rc(int& _rank_mpi, int& _nprow_mpi, int& _npcol_mpi, int& _block_grid, int& _master){

    bool exsist_rc = true;
    FILE* rcfp = NULL;

    //
    // 1. Put .uni10rc file under home directory or same as the binary file.
    // 2. Load the rc file under home directory first.
    //
    
    if(_rank_mpi==0){

      int npnum=1;

      rcfp = fopen("~/.uni10rc", "r");
      if(!rcfp)
        rcfp = fopen(".uni10rc", "r");

      if(!rcfp)
        exsist_rc = false;

      else{

        this->etype = mpi;

        int max_len = 256;
        char buffer[max_len];

        char* pch;
        while(fgets(buffer, max_len, rcfp)){
          
          pch = strtok(buffer, ":");
          pch = strtok(NULL, ":");

          if(strcmp ("NPNUM", buffer)==0)
            npnum = atoi(pch);

          else if(strcmp ("NPROW", buffer)==0)
            _nprow_mpi = atoi(pch);

          else if(strcmp ("NPCOL", buffer)==0)
            _npcol_mpi = atoi(pch);

          else if(strcmp ("BLOCKGRID", buffer)==0)
            _block_grid = atoi(pch);

          else if(strcmp ("MASTERPROC", buffer)==0)
            _master = atoi(pch);
                   
          else
            uni10_error_msg(true, "%s", "Setting the parameters with wrong names.");

        }

        uni10_error_msg(npnum < _nprow_mpi*_npcol_mpi, "%s", " Do not have enough processes available to make a p-by-q process grid.");

      }

      fclose(rcfp);

    }

    MPI_Bcast(&_nprow_mpi , 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&_npcol_mpi , 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&_block_grid, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&_master    , 1, MPI_INT, 0, MPI_COMM_WORLD);

    //printf("row: %d, col: %d, grid: %d, num: %d\n", _nprow_mpi, _npcol_mpi, _block_grid, _rank_mpi);
    MPI_Barrier(MPI_COMM_WORLD);

    return exsist_rc;

  }

  const sysinfo_mpi& uni10_env_mpi::get_info() const{

    return uni10_sys_info;

  }

}
