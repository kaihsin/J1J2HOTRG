#include "uni10/uni10_api/tensor_tools/tensor_tools.h"

namespace uni10{
  
  namespace tensor_tools{

    // Overload for UniTensor<T>::init_para();
    U_para<uni10_double64>* init_para(U_para<uni10_double64>* para, const contain_type style){

      para = init_para_d[style](para);
      return para;

    }

    U_para<uni10_complex128>* init_para(U_para<uni10_complex128>* para, const contain_type style){

      para = init_para_z[style](para);
      return para;

    }

    // Overload for UniTensor<T>::copy_para();
    void copy_para(U_para<uni10_double64>* para, const U_para<uni10_double64>* src_para, const contain_type style){

      copy_para_d[style](para, src_para);

    }

    void copy_para(U_para<uni10_complex128>* para, const U_para<uni10_complex128>* src_para, const contain_type style){

      copy_para_z[style](para, src_para);

    }

    // Overload for UniTensor<T>::init_para();
    void free_para(U_para<uni10_double64>* para, const contain_type style){

      free_para_d[style](para);

    }

    void free_para(U_para<uni10_complex128>* para, const contain_type style){

      free_para_z[style](para);

    }

    // Overload for UniTensor<T>::init();
    void init(U_para<uni10_double64>* para, const contain_type style){

      init_d[style](para);

    }

    void init(U_para<uni10_complex128>* para, const contain_type style){

      init_z[style](para);

    }

    // Overload for UniTensor<T>::initBlocks();
    void initBlocks(U_para<uni10_double64>* para, const contain_type style){

      initBlocks_d[style](para);

    }

    void initBlocks(U_para<uni10_complex128>* para, const contain_type style){

      initBlocks_z[style](para);

    }

    // Overload for UniTensor<T>::setRawElem();
    void setRawElem(U_para<uni10_double64>* para, const uni10_double64* rawElem, const contain_type style){

      setRawElem_d[style](para, rawElem);

    }

    void setRawElem(U_para<uni10_complex128>* para, const uni10_complex128* rawElem, const contain_type style){

      setRawElem_z[style](para, rawElem);

    }

    // Overload for UniTensor<T>::init();
    void putBlock(U_para<uni10_double64>* para, const Qnum& qnum, const Block<uni10_double64>& mat, const contain_type style){

      putBlock_d[style](para, qnum, mat);

    }

    void putBlock(U_para<uni10_complex128>* para, const Qnum& qnum, const Block<uni10_complex128>& mat, const contain_type style){

      putBlock_z[style](para, qnum, mat);

    }

    void set_zero(U_para<uni10_double64>*   para, const contain_type style){

      set_zero_d[style](para);

    }

    void set_zero(U_para<uni10_complex128>* para, const contain_type style){

      set_zero_z[style](para);

    }

    void randomize(U_para<uni10_double64>*   para, const contain_type style){

      randomize_d[style](para);

    }

    void randomize(U_para<uni10_complex128>* para, const contain_type style){

      randomize_z[style](para);

    }

    void transpose(U_para<uni10_double64>* para, const U_para<uni10_double64>*  src_para, const contain_type style){

      transpose_d[style](para, src_para);

    }

    void transpose(U_para<uni10_complex128>* para, const U_para<uni10_complex128>*  src_para, const contain_type style){

      transpose_z[style](para, src_para);

    }

    void dagger(U_para<uni10_double64>* para, const U_para<uni10_double64>*  src_para, const contain_type style){

      dagger_d[style](para, src_para);

    }

    void dagger(U_para<uni10_complex128>* para, const U_para<uni10_complex128>*  src_para, const contain_type style){

      dagger_z[style](para, src_para);

    }

    void conj(U_para<uni10_double64>* para, const U_para<uni10_double64>*  src_para, const contain_type style){

      conj_d[style](para, src_para);

    }

    void conj(U_para<uni10_complex128>* para, const U_para<uni10_complex128>*  src_para, const contain_type style){

      conj_z[style](para, src_para);

    }

    void permute(const U_para<uni10_double64>*   T1_para, const contain_type T1_style, const std::vector<uni10_int32>& rsp_outin, 
        U_para<uni10_double64>*   T2_para, uni10_bool inorder){

      permute_d[T1_style](T1_para, rsp_outin, T2_para, inorder);

    }

    void permute(const U_para<uni10_complex128>* T1_para, const contain_type T1_style, const std::vector<uni10_int32>& rsp_outin, 
        U_para<uni10_complex128>* T2_para, uni10_bool inorder){

      permute_z[T1_style](T1_para, rsp_outin, T2_para, inorder);

    }

    void addGate(U_para<uni10_double64>* para, const std::vector<_Swap>& swaps, const contain_type style){

      addGate_d[style](para, swaps);

    }

    void addGate(U_para<uni10_complex128>* para, const std::vector<_Swap>& swaps, const contain_type style){

      addGate_z[style](para, swaps);

    }

    void traceByRow(U_para<uni10_double64>* Tout_para, const U_para<uni10_double64>* Tin_para, uni10_int32 la, uni10_int32 lb , const contain_type style){

      traceByRow_d[style](Tout_para, Tin_para, la, lb);

    }

    void traceByRow(U_para<uni10_complex128>* Tout_para, const U_para<uni10_complex128>* Tin_para, uni10_int32 la, uni10_int32 lb , const contain_type style){

      traceByRow_z[style](Tout_para, Tin_para, la, lb);

    }

    uni10_double64 tensorAt(U_para<uni10_double64>*   T_para, const uni10_uint64* idxs, const contain_type _style){

      return tensorAt_d[_style](T_para, idxs);

    }

    uni10_complex128 tensorAt(U_para<uni10_complex128>* T_para, const uni10_uint64* idxs, const contain_type _style){

      return tensorAt_z[_style](T_para, idxs);

    }


  };

};
