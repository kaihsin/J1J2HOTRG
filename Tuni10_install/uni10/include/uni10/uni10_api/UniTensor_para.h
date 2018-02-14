#ifndef __UNITENSOR_PARA_H__
#define __UNITENSOR_PARA_H__

#include <stdio.h>

#include <string>
#include <vector>

#include "uni10/uni10_api/Bond.h"
#include "uni10/uni10_api/Block.h"

namespace uni10{

  // Non Symmetry UniTensor
  template<typename uni10_type>
    struct no_sym_para{

      no_sym_para(){};

      // Link to the class UniTensor.
      std::string name;
      std::vector<Bond> bonds;
      std::vector<uni10_int32> labels;
      uni10_int32 RBondNum;                                      //Row bond number
      uni10_uint64 RQdim;
      uni10_uint64 CQdim;
      uni10_uint64 U_elemNum;
      std::map< Qnum, Block<uni10_type> > blocks;
      UELEM(uni10_elem, _package, _type)<uni10_type> U_elem;     // pointer to a real matrix
      uni10_int32 status;                                        // Check initialization, 1 initialized, 3 initialized with label, 5 initialized with elements

    };

  // Has Symmetry and save as block diagonal matrix.
  template<typename uni10_type>
    struct blk_sym_para{

      blk_sym_para(){};

      // Link to the class UniTensor.
      std::string name;
      std::vector<Bond> bonds;
      std::vector<uni10_int32>labels;
      //void packMeta();
      uni10_int32 RBondNum;                                    //Row bond number
      uni10_uint64 RQdim;
      uni10_uint64 CQdim;
      uni10_uint64 U_elemNum;
      std::map<Qnum, Block<uni10_type> > blocks;
      UELEM(uni10_elem, _package, _type)<uni10_type> U_elem;     // pointer to a real matrix
      uni10_int32 status;                                        // Check initialization, 1 initialized, 3 initialized with label, 5 initialized with elements

      // Handle the blk symmetry.
      std::map<uni10_int32, Block<uni10_type>* > RQidx2Blk;    //Qidx to the Block
      std::map<uni10_int32, uni10_uint64> QidxEnc;
      std::map<uni10_int32, uni10_uint64> RQidx2Off;    //the row offset starts from the block origin of a qnum
      std::map<uni10_int32, uni10_uint64> CQidx2Off;    //the col offset starts from the block origin of a qnum
      std::map<uni10_int32, uni10_uint64> RQidx2Dim;
      std::map<uni10_int32, uni10_uint64> CQidx2Dim;

    };

  template<typename uni10_type>
    struct spar_sym_para{

      spar_sym_para(){};
      // Developping 
      std::string name;
      std::vector<Bond> bonds;
      std::map<Qnum, Block<uni10_type> > blocks;
      std::vector<uni10_int32>labels;
      void packMeta();
      uni10_int32 RBondNum;   //Row bond number
      uni10_int32 RQdim;
      uni10_int32 CQdim;
      uni10_uint64 m_elemNum;

    };

  template<typename uni10_type>
    struct U_para{
      U_para(): nsy(NULL), bsy(NULL), ssy(NULL), check_status(-1){};

      void nsy_alloc(){ nsy = new struct no_sym_para<uni10_type>[1]; }

      void free(){
          if(nsy != NULL)
            delete nsy;
          else if(bsy != NULL)
            delete bsy;
          if(ssy != NULL)
            delete ssy;
      }

      no_sym_para<uni10_type>*     nsy;
      blk_sym_para<uni10_type>*    bsy;
      spar_sym_para<uni10_type>*   ssy;

      uni10_int32 check_status;

    };

};

#endif
