#include "uni10/uni10_cusolver_gpu/tools_cusolver_gpu/uni10_tools_cusolver_gpu.h"

namespace uni10{

  void* uni10_elem_alloc(uni10_uint64 memsize, uni10_bool& __ongpu){

    void* ptr = NULL;

    if(RUNTIMETYPE == only_cpu){

      ptr = malloc(memsize);
      __ongpu = false;

    }
    else if(RUNTIMETYPE == hybrid){

      uni10_error_msg(true, "%s", "Developing");

    }
    else if(RUNTIMETYPE == only_gpu){

      checkCudaErrors(cudaMallocManaged(&ptr, memsize));
      __ongpu = true;

    }

    env_variables.used_memsize(memsize);
    uni10_error_msg(ptr==NULL, "%s","Fails in allocating memory.");

    return ptr;
  }

  void* uni10_elem_copy(void* des, const void* src, uni10_uint64 memsize, uni10_bool des_ongpu, uni10_bool src_ongpu){

    // wantgo == 0; Both pointers are on the host.
    // wantgo == 1; The source data pointer is on the device and the destination one is allocated on the host.
    // wantgo == 2; The destination pointer ison the device and the source data one is allocated on the host.
    // wantgo == 3; Both pointers are on the device.
    
    int wantgo = des_ongpu * 2 + src_ongpu; 

    if(wantgo == 0){

      memcpy(des, src, memsize);

    }
    else if(wantgo == 1){

      checkCudaErrors(cudaMemcpy(des, src, memsize, cudaMemcpyDeviceToHost));

    }
    else if(wantgo == 2){

      checkCudaErrors(cudaMemcpy(des, src, memsize, cudaMemcpyHostToDevice));
    }
    else if(wantgo == 3){

      checkCudaErrors(cudaMemcpy(des, src, memsize, cudaMemcpyDeviceToDevice));

    }

    return des;

  }

  void uni10_elem_free(void* ptr, uni10_uint64 memsize, uni10_bool __ongpu){

    if(__ongpu){

      checkCudaErrors(cudaFree(ptr));

    }
    else{

      free(ptr);

    }

    env_variables.used_memsize(-memsize);
    ptr = NULL;

  }

  void uni10_elemBzero(void* ptr, uni10_uint64 memsize, uni10_bool __ongpu){

    if(__ongpu){

      checkCudaErrors(cudaMemset(ptr, 0, memsize));

    }
    else{

      memset(ptr, 0, memsize);

    }

  }

  // For double 
  //
  void uni10_setDiag(uni10_double64* ori_elem, uni10_double64* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n, uni10_bool ori_ongpu, uni10_bool diag_ongpu){

    // wantgo == 0; Both pointers are on the host.
    // wantgo == 1; The pointer allocated for a diagonal matrix is on the device and the other one for rectangular matrix is on the host.
    // wantgo == 2; The pointer allocated for a rectangular matrix is on the host and the other one for a diagonal matrix is on the device.
    // wantgo == 3; Both pointers are on the device.
    
    int wantgo = ori_ongpu * 2 + diag_ongpu; 

    uni10_uint64 min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    if(wantgo == 0){

    for(uni10_uint64 i = 0; i < min; i++)
      ori_elem[i * n + i] = diag_elem[i];

    }
    else{

      uni10_error_msg(true, "%s", "Developing");

    }

  }

  void uni10_getDiag(uni10_double64* ori_elem, uni10_double64* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n, uni10_bool ori_ongpu, uni10_bool diag_ongpu){

    // wantgo == 0; Both pointers are on the host.
    // wantgo == 1; The pointer allocated for a diagonal matrix is on the device and the other one for rectangular matrix is on the host.
    // wantgo == 2; The pointer allocated for a rectangular matrix is on the host and the other one for a diagonal matrix is on the device.
    // wantgo == 3; Both pointers are on the device.
    
    int wantgo = ori_ongpu * 2 + diag_ongpu; 

    uni10_uint64 min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;
    
    if(wantgo == 0){

    for(uni10_uint64 i = 0; i < min; i++)
      diag_elem[i] = ori_elem[i * n + i];

    }
    else{

      uni10_error_msg(true, "%s", "Developing");

    }

  }

  void uni10_getUpTri(uni10_double64* ori_elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n, uni10_bool ori_ongpu, uni10_bool tri_ongpu){

    // wantgo == 0; Both pointers are on the host.
    // wantgo == 1; The pointer allocated for a upper triangular matrix is on the device and the other one for rectangular matrix is on the host.
    // wantgo == 2; The pointer allocated for a rectangular matrix is on the host and the other one for a upper triangular matrix is on the device.
    // wantgo == 3; Both pointers are on the device.
    
    int wantgo = ori_ongpu * 2 + tri_ongpu; 

    uni10_uint64 min = m < n ? m : n;

    if(wantgo == 0){

      for(uni10_uint64 i = 0; i < min; i++)
        uni10_elem_copy(tri_elem + i*min+i, ori_elem + i*n + (n-min)+i, (min-i)*sizeof(uni10_double64), false, false);

    }
    else{

      uni10_error_msg(true, "%s", "Developing");

    }

  }

  void uni10_getDnTri(uni10_double64* ori_elem, uni10_double64* tri_elem, uni10_uint64 m, uni10_uint64 n, uni10_bool ori_ongpu, uni10_bool tri_ongpu){

    // wantgo == 0; Both pointers are on the host.
    // wantgo == 1; The pointer allocated for a lower triangular matrix is on the device and the other one for rectangular matrix is on the host.
    // wantgo == 2; The pointer allocated for a rectangular matrix is on the host and the other one for a lower triangular matrix is on the device.
    // wantgo == 3; Both pointers are on the device.

    int wantgo = ori_ongpu * 2 + tri_ongpu; 

    uni10_uint64 min = m < n ? m : n;

    if(wantgo == 0){

    for(uni10_uint64 i = 0; i < min; i++)
      uni10_elem_copy(tri_elem + i*min, ori_elem + n*(m-min)+i*n, (i+1)*sizeof(uni10_double64), false, false);

    }
    else{

      uni10_error_msg(true, "%s", "Developing");

    }

  }

  void uni10_print_elem_i(const uni10_double64& elem_i){

    fprintf(stdout, " %8.4f", elem_i);

  }

  uni10_double64 UNI10_REAL( uni10_double64 elem_i ){

    return elem_i;

  }

  uni10_double64 UNI10_IMAG( uni10_double64 elem_i ){

    return 0.;
   
  }

  // For complex 
  //
  void uni10_setDiag(uni10_complex128* ori_elem, uni10_complex128* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n, uni10_bool ori_ongpu, uni10_bool diag_ongpu){

    // wantgo == 0; Both pointers are on the host.
    // wantgo == 1; The pointer allocated for a diagonal matrix is on the device and the other one for rectangular matrix is on the host.
    // wantgo == 2; The pointer allocated for a rectangular matrix is on the host and the other one for a diagonal matrix is on the device.
    // wantgo == 3; Both pointers are on the device.
    
    int wantgo = ori_ongpu * 2 + diag_ongpu; 

    uni10_uint64 min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    if(wantgo == 0){

    for(uni10_uint64 i = 0; i < min; i++)
      ori_elem[i * n + i] = diag_elem[i];

    }
    else{

      uni10_error_msg(true, "%s", "Developing");

    }

  }

  void uni10_getDiag(uni10_complex128* ori_elem, uni10_complex128* diag_elem, uni10_uint64 m, uni10_uint64 n, uni10_uint64 diag_n, uni10_bool ori_ongpu, uni10_bool diag_ongpu){

    // wantgo == 0; Both pointers are on the host.
    // wantgo == 1; The pointer allocated for a diagonal matrix is on the device and the other one for rectangular matrix is on the host.
    // wantgo == 2; The pointer allocated for a rectangular matrix is on the host and the other one for a diagonal matrix is on the device.
    // wantgo == 3; Both pointers are on the device.
    
    int wantgo = ori_ongpu * 2 + diag_ongpu; 

    uni10_uint64 min = m < n ? m : n;

    min = min < diag_n ? min : diag_n;

    if(wantgo == 0){

    for(uni10_uint64 i = 0; i < min; i++)
      diag_elem[i] = ori_elem[i * n + i];

    }
    else{

      uni10_error_msg(true, "%s", "Developing");

    }

  }

  void uni10_getUpTri(uni10_complex128* ori_elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n, uni10_bool ori_ongpu, uni10_bool tri_ongpu){

    // wantgo == 0; Both pointers are on the host.
    // wantgo == 1; The pointer allocated for a upper triangular matrix is on the device and the other one for rectangular matrix is on the host.
    // wantgo == 2; The pointer allocated for a rectangular matrix is on the host and the other one for a upper triangular matrix is on the device.
    // wantgo == 3; Both pointers are on the device.
    
    int wantgo = ori_ongpu * 2 + tri_ongpu; 

    uni10_uint64 min = m < n ? m : n;

    if(wantgo == 0){

      for(uni10_uint64 i = 0; i < min; i++)
        uni10_elem_copy(tri_elem + i*min+i, ori_elem + i*n+(n-min)+i, (min-i)*sizeof(uni10_complex128), false, false);

    }
    else{

      uni10_error_msg(true, "%s", "Developing");

    }

  }

  void uni10_getDnTri(uni10_complex128* ori_elem, uni10_complex128* tri_elem, uni10_uint64 m, uni10_uint64 n, uni10_bool ori_ongpu, uni10_bool tri_ongpu){

    // wantgo == 0; Both pointers are on the host.
    // wantgo == 1; The pointer allocated for a lower triangular matrix is on the device and the other one for rectangular matrix is on the host.
    // wantgo == 2; The pointer allocated for a rectangular matrix is on the host and the other one for a lower triangular matrix is on the device.
    // wantgo == 3; Both pointers are on the device.

    int wantgo = ori_ongpu * 2 + tri_ongpu; 

    uni10_uint64 min = m < n ? m : n;

    if(wantgo == 0){

    for(uni10_uint64 i = 0; i < min; i++)
      uni10_elem_copy(tri_elem + i*min, ori_elem + n*(m-min)+i*n, (i+1)*sizeof(uni10_complex128), false, false);

    }
    else{

      uni10_error_msg(true, "%s", "Developing");

    }

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
  void uni10_elem_cast(uni10_complex128* des, const uni10_double64* src, uni10_uint64 N, bool des_ongpu, bool src_ongpu){

    uni10_error_msg(true, "%s", "Developing");

  }

  void uni10_elem_cast(uni10_double64* des, const uni10_complex128* src, uni10_uint64 N, bool des_ongpu, bool src_ongpu){

    uni10_error_msg(true, "%s", "Developing");

  }

  void shrinkWithoutFree(uni10_uint64 memsize){

    env_variables.used_memsize(-memsize);

  }

} /* namespace uni10 */
