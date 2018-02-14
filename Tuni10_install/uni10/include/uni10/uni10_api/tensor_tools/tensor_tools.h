#ifndef __UNI10_TENSOR_TOOLS_H__
#define __UNI10_TENSOR_TOOLS_H__

#include <stdio.h>

#include "uni10/uni10_api/tensor_tools/nsy_tensor_tools.h"
#include "uni10/uni10_api/tensor_tools/bsy_tensor_tools.h"

namespace uni10{

  namespace tensor_tools{
    // The function tools to initialize UniTensor.
    static U_para<uni10_double64>*  (*init_para_d[])(U_para<uni10_double64>*  ) = {init_para_nsy, init_para_bsy};
    static U_para<uni10_complex128>*(*init_para_z[])(U_para<uni10_complex128>*) = {init_para_nsy, init_para_bsy};

    // The function tools to copy parameters.
    static void (*copy_para_d[])(U_para<uni10_double64>*  , const U_para<uni10_double64>*   ) = {copy_para_nsy, copy_para_bsy};
    static void (*copy_para_z[])(U_para<uni10_complex128>*, const U_para<uni10_complex128>* ) = {copy_para_nsy, copy_para_bsy};

    // The function tools to initialize UniTensor.
    static void (*free_para_d[])(U_para<uni10_double64>*  ) = {free_para_nsy, free_para_bsy};
    static void (*free_para_z[])(U_para<uni10_complex128>*) = {free_para_nsy, free_para_bsy};

    static void (*init_d[])(U_para<uni10_double64>*  ) = {init_nsy, init_bsy};
    static void (*init_z[])(U_para<uni10_complex128>*) = {init_nsy, init_bsy};

    static void (*initBlocks_d[])(U_para<uni10_double64>*  ) = {initBlocks_nsy, initBlocks_bsy};
    static void (*initBlocks_z[])(U_para<uni10_complex128>*) = {initBlocks_nsy, initBlocks_bsy};

    static void (*setRawElem_d[])(U_para<uni10_double64>*   T1_para, const uni10_double64*   rawElem) = {setRawElem_nsy, setRawElem_bsy};
    static void (*setRawElem_z[])(U_para<uni10_complex128>* T1_para, const uni10_complex128* rawElem) = {setRawElem_nsy, setRawElem_bsy};

    static void (*putBlock_d[])(U_para<uni10_double64>*  ,const Qnum& qnum, const Block<uni10_double64>& mat  ) = {putBlock_nsy, putBlock_bsy};
    static void (*putBlock_z[])(U_para<uni10_complex128>*,const Qnum& qnum, const Block<uni10_complex128>& mat) = {putBlock_nsy, putBlock_bsy};

    static void (*set_zero_d[])(U_para<uni10_double64>*  ) = {set_zero_nsy, set_zero_bsy};
    static void (*set_zero_z[])(U_para<uni10_complex128>*) = {set_zero_nsy, set_zero_bsy};

    static void (*randomize_d[])(U_para<uni10_double64>*  ) = {randomize_nsy, randomize_bsy};
    static void (*randomize_z[])(U_para<uni10_complex128>*) = {randomize_nsy, randomize_bsy};

    static void (*transpose_d[])(U_para<uni10_double64>*  , const U_para<uni10_double64>*  ) = {transpose_nsy, transpose_bsy};
    static void (*transpose_z[])(U_para<uni10_complex128>*, const U_para<uni10_complex128>*) = {transpose_nsy, transpose_bsy};

    static void (*dagger_d[])(U_para<uni10_double64>*  , const U_para<uni10_double64>*  ) = {dagger_nsy, dagger_bsy};
    static void (*dagger_z[])(U_para<uni10_complex128>*, const U_para<uni10_complex128>*) = {dagger_nsy, dagger_bsy};

    static void (*conj_d[])(U_para<uni10_double64>*  , const U_para<uni10_double64>*  ) = {conj_nsy, conj_bsy};
    static void (*conj_z[])(U_para<uni10_complex128>*, const U_para<uni10_complex128>*) = {conj_nsy, conj_bsy};

    static void (*permute_d[])(const U_para<uni10_double64>*   T1_para, const std::vector<uni10_int32>& rsp_outin, 
        U_para<uni10_double64>*   T2_para, uni10_bool inorder) = {permute_nsy, permute_bsy};
    static void (*permute_z[])(const U_para<uni10_complex128>* T1_para, const std::vector<uni10_int32>& rsp_outin,
        U_para<uni10_complex128>* T2_para, uni10_bool inorder) = {permute_nsy, permute_bsy};

    static void (*addGate_d[])(U_para<uni10_double64>* T1_para  , const std::vector<_Swap>& swaps) = {addGate_nsy, addGate_bsy};
    static void (*addGate_z[])(U_para<uni10_complex128>* T1_para, const std::vector<_Swap>& swaps) = {addGate_nsy, addGate_bsy};

    static void (*traceByRow_d[])(U_para<uni10_double64>* Tout_para  , const U_para<uni10_double64>* Tin_para  , uni10_int32 la, uni10_int32 lb ) = {traceByRow_nsy, traceByRow_bsy};
    static void (*traceByRow_z[])(U_para<uni10_complex128>* Tout_para, const U_para<uni10_complex128>* Tin_para, uni10_int32 la, uni10_int32 lb ) = {traceByRow_nsy, traceByRow_bsy};

    static uni10_double64   (*tensorAt_d[])(U_para<uni10_double64>*   T_para, const uni10_uint64* idxs) = {tensorAt_nsy, tensorAt_bsy};
    static uni10_complex128 (*tensorAt_z[])(U_para<uni10_complex128>* T_para, const uni10_uint64* idxs) = {tensorAt_nsy, tensorAt_bsy};



    // Function overload for initialize U_paras;
    U_para<uni10_double64>*   init_para(U_para<uni10_double64>*   const para, contain_type _style);
    U_para<uni10_complex128>* init_para(U_para<uni10_complex128>* const para, contain_type _style);

    // Function overload for free parameters;
    void copy_para(U_para<uni10_double64>*   para, const U_para<uni10_double64>*   src_para, const contain_type _style);
    void copy_para(U_para<uni10_complex128>* para, const U_para<uni10_complex128>* src_para, const contain_type _style);

    // Function overload for free parameters;
    void free_para(U_para<uni10_double64>*   para, const contain_type _style);
    void free_para(U_para<uni10_complex128>* para, const contain_type _style);

    // Function overload for UniTensor<T>::init();
    void init(U_para<uni10_double64>*   para, const contain_type _style);
    void init(U_para<uni10_complex128>* para, const contain_type _style);
   
    // Function overload for UniTensor<T>::initBlocks();
    void initBlocks(U_para<uni10_double64>*   para, const contain_type _style);
    void initBlocks(U_para<uni10_complex128>* para, const contain_type _style);

    void setRawElem(U_para<uni10_double64>*   para, const uni10_double64*   rawElem, const contain_type _style);
    void setRawElem(U_para<uni10_complex128>* para, const uni10_complex128* rawElem, const contain_type _style);

    void putBlock(U_para<uni10_double64>*   para, const Qnum& qnum, const Block<uni10_double64>& mat  , const contain_type _style);
    void putBlock(U_para<uni10_complex128>* para, const Qnum& qnum, const Block<uni10_complex128>& mat, const contain_type _style);

    void setElem(U_para<uni10_double64>*   para, const Qnum& qnum, const Block<uni10_double64>& mat  , const contain_type _style);
    void setElem(U_para<uni10_complex128>* para, const Qnum& qnum, const Block<uni10_complex128>& mat, const contain_type _style);

    void set_zero(U_para<uni10_double64>*   para, const contain_type _style);
    void set_zero(U_para<uni10_complex128>* para, const contain_type _style);

    void randomize(U_para<uni10_double64>*   para, const contain_type _style);
    void randomize(U_para<uni10_complex128>* para, const contain_type _style);

    void randomize(U_para<uni10_double64>*   para, const contain_type _style);
    void randomize(U_para<uni10_complex128>* para, const contain_type _style);

    void transpose(U_para<uni10_double64>*   para,const U_para<uni10_double64>*   src_para, const contain_type style);
    void transpose(U_para<uni10_complex128>* para,const U_para<uni10_complex128>* src_para, const contain_type style);

    void dagger(U_para<uni10_double64>*   para,const U_para<uni10_double64>*   src_para, const contain_type style);
    void dagger(U_para<uni10_complex128>* para,const U_para<uni10_complex128>* src_para, const contain_type style);

    void conj(U_para<uni10_double64>*   para,const U_para<uni10_double64>*   src_para, const contain_type style);
    void conj(U_para<uni10_complex128>* para,const U_para<uni10_complex128>* src_para, const contain_type style);

    void permute(const U_para<uni10_double64>*   T1_para, const contain_type T1_style, const std::vector<uni10_int32>& rsp_outin, 
        U_para<uni10_double64>*   T2_para, uni10_bool inorder);
    void permute(const U_para<uni10_complex128>* T1_para, const contain_type T1_style, const std::vector<uni10_int32>& rsp_outin,
        U_para<uni10_complex128>* T2_para, uni10_bool inorder);

    void addGate(U_para<uni10_double64>*   T1_para, const std::vector<_Swap>& swaps, const contain_type _style);
    void addGate(U_para<uni10_complex128>* T1_para, const std::vector<_Swap>& swaps, const contain_type _style);

    void traceByRow(U_para<uni10_double64>*   Tout_para, const U_para<uni10_double64>*   Tin_para, uni10_int32 la, uni10_int32 lb , const contain_type _style);
    void traceByRow(U_para<uni10_complex128>* Tout_para, const U_para<uni10_complex128>* Tin_para, uni10_int32 la, uni10_int32 lb , const contain_type _style);

    uni10_double64   tensorAt(U_para<uni10_double64>*   T_para, const uni10_uint64* idxs, const contain_type _style);
    uni10_complex128 tensorAt(U_para<uni10_complex128>* T_para, const uni10_uint64* idxs, const contain_type _style);

  };

};

#endif
