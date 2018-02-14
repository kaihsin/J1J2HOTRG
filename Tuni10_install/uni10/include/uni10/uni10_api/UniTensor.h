/****************************************************************************
 *  @file UniTensor.h
 *  @license
 *   Universal Tensor Network Library
 *   Copyright (c) 2013-2014
 *   National Taiwan University
 *   National Tsing-Hua University
 *
 *   This file is part of Uni10, the Universal Tensor Network Library.
 *
 *   Uni10 is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Uni10 is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
 *  @endlicense
 *  @brief Header file for UniTensor class
 *  @author Yun-Da Hsieh
 *  @author Yun-Hsuan Chou
 *  @date 2014-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/
#ifndef __UNI10_UNITENSOR_H__
#define __UNI10_UNITENSOR_H__

#include <string>

#include "uni10/uni10_api/Matrix.h"
#include "uni10/uni10_api/UniTensor_para.h"
#include "uni10/uni10_api/network_tools/network_tools.h"

/// @brief Uni10 - the Universal Tensor %Network Library
namespace uni10{

  enum contain_type{
    no_sym    = 0,
    blk_sym   = 1,
    spar_sym  = 2 
  };

  struct nsy_paras;

  template<typename uni10_type>
    class UniTensor;

  template<typename uni10_type>
    class Node;

  template<typename uni10_type>
    class Network;

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, const UniTensor<uni10_type>& _b);

  template<typename uni10_type>
    UniTensor<uni10_type> operator*(const UniTensor<uni10_type>& Ta, uni10_double64 a);

  template<typename uni10_type>
    UniTensor<uni10_complex128> operator*(const UniTensor<uni10_type>& Ta, uni10_complex128 a);

  template<typename uni10_type>
    UniTensor<uni10_type> operator*(uni10_double64 a, const UniTensor<uni10_type>& Ta);

  template<typename uni10_type>
    UniTensor<uni10_complex128> operator*(uni10_complex128 a, const UniTensor<uni10_type>& Ta);

  // UniTensor element-wise addition.
  template<typename T>
    UniTensor<T> operator+(const UniTensor<uni10_double64>& t1, const UniTensor<T>& t2);

  template<typename T>
    UniTensor<uni10_complex128> operator+(const UniTensor<uni10_complex128>& t1, const UniTensor<T>& t2);

  // UniTensor element-wise subtract.
  template<typename T>
    UniTensor<T> operator-(const UniTensor<uni10_double64>& t1, const UniTensor<T>& t2);

  template<typename T>
    UniTensor<uni10_complex128> operator-(const UniTensor<uni10_complex128>& t1, const UniTensor<T>& t2);

  UniTensor<uni10_complex128> operator-(const UniTensor<uni10_double64>& t1, const UniTensor<uni10_complex128>& t2);

  UniTensor<uni10_complex128>& operator+=(UniTensor<uni10_complex128>& _m1, const UniTensor<uni10_double64>& _m2);

  UniTensor<uni10_complex128>& operator-=(UniTensor<uni10_complex128>& _m1, const UniTensor<uni10_double64>& _m2);

  UniTensor<uni10_complex128>& operator*=(UniTensor<uni10_complex128>& _m1, uni10_complex128 a);

  ///@class UniTensor
  ///@brief The UniTensor class defines the symmetric tensors
  ///
  /// A UniTensor consists of Bond's carrying quantum numbers Qnum's. The tensor elements are organized as
  /// quantum number blocks. The Qnum's on the Bonds defines the size of the Qnum blocks and the rank of
  /// UniTensor is defined by the number of Bond's.\par
  /// Each Bond carries a label. Labels are used to manipulate tensors, such as in permute, partialTrace and
  /// contraction. \par
  /// Operations on tensor elements is pefromed through  getBlock and putBlock functions to take out/put in
  /// block elements out as a Matrix.
  /// @see Qnum, Bond, Matrix< uni10_type >

  template <typename uni10_type>
    class UniTensor {

      private:

        contain_type style;
        struct U_para<uni10_type>* paras;

        void init_para();
        void meta_link();
        void copy_para(U_para<uni10_type>* src_para);
        void init();
        void check_status();
        void initBlocks();
        void free_para();
        void init_paras_null();
        contain_type check_bonds(const std::vector<Bond>& _bonds)const;

        // General variables.
        std::string* name;
        std::vector<Bond>* bonds;
        std::vector<uni10_int>* labels;
        uni10_int*  RBondNum;       //Row bond number
        uni10_uint64* RQdim;
        uni10_uint64* CQdim;
        uni10_uint64* U_elemNum;
        std::map< Qnum, Block<uni10_type> >* blocks;
        UELEM(uni10_elem, _package, _type)<uni10_type>* U_elem;     // pointer to a real matrix
        uni10_int*  status;     //Check initialization, 1 initialized, 3 initialized with label, 5 initialized with elements

        static uni10_int COUNTER;
        static uni10_uint64 ELEMNUM;
        static uni10_uint64 MAXELEMNUM;
        static uni10_uint64 MAXELEMTEN;               //Max number of element of a tensor

        static const uni10_int HAVEBOND = 1;        /**< A flag for initialization */
        static const uni10_int HAVEELEM = 2;        /**< A flag for having element assigned */

      public:

        static uni10_int GET_HAVEBOND(){return HAVEBOND;}; 
        static uni10_int GET_HAVEELEM(){return HAVEELEM;};      /**< A flag for having element assigned */
        static uni10_int GET_COUNTER(){return COUNTER;}
        static uni10_uint64 GET_ELEMNUM(){return ELEMNUM;}
        static uni10_uint64 GET_MAXELEMNUM(){return MAXELEMNUM;};
        static uni10_uint64 GET_MAXELEMTEN(){return MAXELEMTEN;};

        ///
        /// @brief Default Constructor
        ///
        explicit UniTensor();
        /// @brief Create a rank-0 UniTensor with value \c val
        ///
        /// @param val Value of the scalar
        explicit UniTensor(uni10_type val);

        /// @brief Create a UniTensor from a list of Bond's
        /// @param _bonds List of bonds
        /// @param _name Name of the tensor, defaults to ""
        explicit UniTensor(const std::vector<Bond>& _bonds, const std::string& _name = "");

        /// @brief Create a UniTensor from a list of Bond's and assign labels
        /// @param _bonds List of bonds
        /// @param labels Labels for \c _bonds
        /// @param _name Name of the tensor, defaults to ""
        explicit UniTensor(const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");

        /// @brief Create a UniTensor from a list of Bond's and assign labels
        /// @param _bonds List of bonds
        /// @param labels Labels for \c _bonds
        /// @param _name Name of the tensor, defaults to ""
        explicit UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");

        //explicit UniTensor(const std::string& fname);
       
        /// @brief Create a UniTensor from a Block< uni10_type >
        explicit UniTensor(const Block<uni10_type>& UniT);

        /// @brief Copy constructor
        UniTensor(const UniTensor& UniT);

        template <typename U>
          UniTensor(const UniTensor<U>& UniT);

        /// @brief Copy constructor
        UniTensor(const std::string& fname);

        /// @brief Destructor
        ///
        /// Destroys the UniTensor and freeing all allocated memory
        ~UniTensor();

        void save(const std::string& fname) const;

        void load(const std::string& fname);

        void set_zeros();

        void set_zeros(const Qnum& qnum);

        void identity();

        void identity(const Qnum& qnum);

        void randomize(char UorN='U', uni10_double64 dn_mu=-1, uni10_double64 up_var=1, uni10_int64 seed=uni10_clock);

        void randomize(const Qnum& qnum, char UorN='U', uni10_double64 dn_mu=-1, uni10_double64 up_var=1, uni10_int64 seed=uni10_clock );

        void orthoRand(char UorN='U', uni10_double64 dn_mu=-1, uni10_double64 up_var=1, uni10_int64 seed=uni10_clock);

        void orthoRand(const Qnum& qnum, char UorN='U', uni10_double64 dn_mu=-1, uni10_double64 up_var=1, uni10_int64 seed=uni10_clock );

        /// @brief Assign elements to a block
        ///
        /// Assigns elements of the  matrix \c mat to the  block of quantum number \c qnum, replacing the origin
        /// elements. \par
        /// If \c mat is diagonal,  all the off-diagonal elements are set to zero.
        /// @param mat The matrix elements to be assigned
        void putBlock(const Block<uni10_type>& mat);

        /// @brief Assign elements to a block
        ///
        /// Assigns elements of the  matrix \c mat to the  Qnum(0) block, for non-symmetry tensors.
        ///
        /// If \c mat is diagonal,  all the off-diagonal elements are set to zero.
        /// @param qnum quantum number of the block
        /// @param mat The matrix elements to be assigned
        void putBlock(const Qnum& qnum, const Block<uni10_type>& mat);
        
        /// @brief Access labels
        ///
        /// Returns the labels of the bonds in UniTensor< uni10_type >
        ///
        /// @return List of labels
        std::vector<int> label()const{return *labels;};

        /// @brief Access label
        ///
        /// Access the label of Bond \c idx
        /// @param idx Bond index
        /// @return Label of Bond \c idx
        uni10_int label(size_t idx)const{return (*labels)[idx];};

        /// @brief Access name
        ///
        /// Return the name of the UniTensor< uni10_type >.
        std::string getName() const{return *name;};

        /// @brief Assign name
        ///
        /// Assigns name to the UniTensor< uni10_type >.
        /// @param name Name to be assigned
        void setName(const std::string& _name);

        ///
        const std::map<Qnum, Block<uni10_type> >& const_getBlocks()const{return *blocks;};

        /// @brief Access the number of bonds
        ///
        /// @return Number of bonds
        uni10_uint64 bondNum()const{return bonds->size();};

        /// @brief Access the number of incoming bonds
        ///
        /// @return Number of incoming bonds
        uni10_uint64 inBondNum()const{return *RBondNum;};

        /// @brief Access bonds
        ///
        /// Returns the bonds in UniTensor
        /// @return List of bonds
        std::vector<Bond> bond()const{return *bonds;};

        /// @brief Access bond
        ///
        /// Returns the bond at the position \c idx
        /// @param idx Position of the bond being retrieved
        /// @return A bond
        Bond bond(uni10_uint64 idx)const{return (*bonds)[idx];};

        /// @brief Access the number of elements
        ///
        /// Returns the number of total elements of the blocks.
        /// @return  Number of elements
        uni10_uint64 elemNum()const{return (*U_elemNum); };

        uni10_uint64 typeID()const{return (U_elem->__uni10_typeid); };

        /// @brief Access the number of blocks
        ///
        /// Returns the number of blocks
        /// @return The number of blocks
        uni10_uint64 blockNum()const{return blocks->size();};

        ///
        uni10_type* getElem()const{return U_elem->__elem; }

        Matrix<uni10_type> getRawElem()const;

        /// @brief Get elements in a block
        ///
        /// Returns elements of Qnum(0) Block.
        /// @return The Qnum(0) Block
        const Block<uni10_type>& const_getBlock()const;

        /// @brief Get elements in a block
        ///
        /// Returns elements of \c qnum block.
        /// @return The \c qnum  Block
        const Block<uni10_type>& const_getBlock(const Qnum& qnum)const;

        ///
        void setLabel(const uni10_int newLabel, const uni10_uint64 idx);

        /// @brief Assign labels to bonds in UniTensor
        ///
        /// Assigns the labels \c newLabels to each bond of  UniTensor, replacing the origin labels on the bonds.
        /// @param newLabels Array of labels
        void setLabel(const std::vector<uni10_int>& newLabels);

        /// @overload
        void setLabel(uni10_int* newLabels);

        /// @brief Access block quantum numbers
        ///
        /// Returns the quantum numbers for all blocks in UniTensor.
        /// The return array of quantum numbers is in the ascending order defined in Qnum.
        /// @return Array of Qnum's
        std::vector<Qnum> blockQnum()const;

        /// @brief Access block quantum number
        ///
        /// Returns the quantum number for block \c idx in UniTensor.
        /// Blocks are orderd in the ascending order of Qnum
        /// @param idx Block index
        /// @return Quantum number of block \c idx
        Qnum blockQnum(uni10_uint64 idx)const;

        std::map< Qnum, Matrix<uni10_type> > getBlocks()const;

        Matrix<uni10_type> getBlock(bool diag = false)const;

        Matrix<uni10_type> getBlock(const Qnum& qnum, bool diag = false)const;

        /// @brief Assign raw elements
        ///
        /// Assigns raw elements in \c blk to UniTensor.
        ///
        /// This function will reorganize the raw elements into the block-diagonal form.
        /// @param blk  Block of raw elements
        void setRawElem(const std::vector<uni10_type>& rawElem);

        /// @overload
        void setRawElem(const uni10_type* rawElem);

        /// @overload
        void setRawElem(const Block<uni10_type>& blk);

        void setElem(const uni10_type* _elem);

        void setElem(const std::vector<uni10_type>& _elem);

        UniTensor& assign(const std::vector<Bond>& _bond);

        UniTensor& combineBond(const std::vector<uni10_int>& combined_labels);

        UniTensor& combineBond(uni10_int* combined_labels, uni10_int boundNum);

        uni10_type at(const std::vector<uni10_uint64>& idxs) const;

        void addGate(const std::vector<_Swap>& swaps);

        std::vector<_Swap> exSwap(const UniTensor<uni10_double64>& Tb)const;

        std::vector<_Swap> exSwap(const UniTensor<uni10_complex128>& Tb)const;

        /// @brief Print the diagrammatic representation of UniTensor< uni10_type >
        ///
        /// Prints out the diagrammatic representation of UniTensor \c uT as (for example):
        /// @code
        ///**************** Demo ****************
        ///     ____________
        ///    |            |
        ///0___|3          3|___2
        ///    |            |
        ///1___|3          3|___3
        ///    |            |
        ///    |____________|
        ///
        ///**************************************
        /// @endcode
        void printDiagram()const;

        std::string printRawElem(bool print=true)const;

        static std::string profile(bool print = true);

        void clear();

        UniTensor& operator=(UniTensor const& UniT);

        template<typename U>
          UniTensor& operator=(UniTensor<U> const& UniT);

        UniTensor<uni10_type>& operator+=( const UniTensor<uni10_type>& Tb );

        UniTensor<uni10_type>& operator-=( const UniTensor<uni10_type>& Tb );

        UniTensor& operator*= (uni10_double64 a);

        //friend UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb);
        uni10_type operator[](uni10_uint64 idx)const{
          uni10_error_msg(!(idx < (*U_elemNum)), "Index exceeds the number of elements( %ld ).", *U_elemNum);
          return U_elem->__elem[idx];
        }

        /// @brief Print out UniTensor
        ///
        /// Prints out a UniTensor \c uT as(for example):
        /// @code
        ///**************** Demo ****************
        ///     ____________
        ///    |            |
        ///0___|3          3|___2
        ///    |            |
        ///1___|3          3|___3
        ///    |            |
        ///    |____________|
        ///
        ///================BONDS===============
        ///IN : (U1 = 1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|1, (U1 = -1, P = 0, 0)|1, Dim = 3
        ///IN : (U1 = 1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|1, (U1 = -1, P = 0, 0)|1, Dim = 3
        ///OUT: (U1 = 1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|1, (U1 = -1, P = 0, 0)|1, Dim = 3
        ///OUT: (U1 = 1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|1, (U1 = -1, P = 0, 0)|1, Dim = 3
        ///
        ///===============BLOCKS===============
        ///--- (U1 = -2, P = 0, 0): 1 x 1 = 1
        ///
        ///0.840
        ///
        ///--- (U1 = -1, P = 0, 0): 2 x 2 = 4
        ///
        ///0.394  0.783
        ///
        ///0.798  0.912
        ///
        ///--- (U1 = 0, P = 0, 0): 3 x 3 = 9
        ///
        ///0.198  0.335  0.768
        ///
        ///0.278  0.554  0.477
        ///
        ///0.629  0.365  0.513
        ///
        ///--- (U1 = 1, P = 0, 0): 2 x 2 = 4
        ///
        ///0.952  0.916
        ///
        ///0.636  0.717
        ///
        ///--- (U1 = 2, P = 0, 0): 1 x 1 = 1
        ///
        ///0.142
        ///
        ///Total elemNum: 19
        ///***************** END ****************
        /// @endcode
        ///  In the above example, \c uT has four bonds with default labels [0, 1, 2, 3]. The bonds 0 and 1 are
        /// incoming bonds, and  2, 3 are out-going bonds. Each bond has three states
        /// corresponding to three U1 quantum number [-1, 0, 1]. The block elements of the
        /// tensor are als shown. There are five blocks of various <tt>U1= [-2, -1, 0, 1, 2]</tt> and
        /// various sizes. The total element number is 19.
        friend std::ostream& operator<< <>(std::ostream& os, const UniTensor& _b);  // --> uni10_elem().print_elem()

        template<typename _uni10_type>
          friend UniTensor<_uni10_type> operator*(const UniTensor<_uni10_type>& Ta, uni10_double64 a);

        template<typename _uni10_type>
          friend UniTensor<uni10_complex128> operator*(const UniTensor<_uni10_type>& Ta, uni10_complex128 a);

        template<typename _uni10_type>
          friend UniTensor<_uni10_type> operator*(uni10_double64 a, const UniTensor<_uni10_type>& Ta);

        template<typename _uni10_type>
          friend UniTensor<uni10_complex128> operator*(uni10_complex128 a, const UniTensor<_uni10_type>& Ta);

        // UniTensor element-wise addition.
        template<typename _T>
          friend UniTensor<_T> operator+(const UniTensor<uni10_double64>& t1, const UniTensor<_T>& t2);

        template<typename _T>
          friend UniTensor<uni10_complex128> operator+(const UniTensor<uni10_complex128>& t1, const UniTensor<_T>& t2);

        // UniTensor element-wise subtract.
        template<typename _T>
          friend UniTensor<_T> operator-(const UniTensor<uni10_double64>& t1, const UniTensor<_T>& t2);

        template<typename _T>
          friend UniTensor<uni10_complex128> operator-(const UniTensor<uni10_complex128>& t1, const UniTensor<_T>& t2);

        friend UniTensor<uni10_complex128>& operator+=(UniTensor<uni10_complex128>& _m1, const UniTensor<uni10_double64>& _m2);

        friend UniTensor<uni10_complex128>& operator-=(UniTensor<uni10_complex128>& _m1, const UniTensor<uni10_double64>& _m2);

        friend UniTensor<uni10_complex128>& operator*=(UniTensor<uni10_complex128>& _m1, uni10_complex128 a);
        // The prototypes of high-rank linear algebras.
        template<typename _uni10_type>
          friend void permute( UniTensor<_uni10_type>& Tout, const UniTensor<_uni10_type>& T, const std::vector<uni10_int>& newLabels, uni10_int inBondNum, UNI10_INPLACE on);

        template<typename _uni10_type>
          friend void permute( UniTensor<_uni10_type>& T, const std::vector<uni10_int>& newLabels, uni10_int inBondNum, UNI10_INPLACE on);

        template<typename _To, typename _T, typename _U>
          friend void contract( UniTensor<_To>& Tc, UniTensor<_T>& Ta, UniTensor<_U>& Tb, bool fast, UNI10_INPLACE on);

        template<typename _uni10_type>
          friend void transpose( UniTensor<_uni10_type>& Tout, const UniTensor<_uni10_type>& Tin, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void dagger( UniTensor<_uni10_type>& Tout, const UniTensor<_uni10_type>& Tin, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void conj( UniTensor<_uni10_type>& Tout, const UniTensor<_uni10_type>& Tin, UNI10_INPLACE on );

        template<typename _uni10_type>
          friend void partialTrace(UniTensor<_uni10_type>& Tout, const UniTensor<_uni10_type>& Tin, uni10_int la, uni10_int lb, UNI10_INPLACE on);

        template<typename _uni10_type>
          friend UniTensor<_uni10_type> permute( const UniTensor<_uni10_type>& T, const std::vector<uni10_int>& newLabels, uni10_int inBondNum);

        template<typename _uni10_type>
          friend UniTensor<_uni10_type> permute( const UniTensor<_uni10_type>& T, uni10_int* newLabels, uni10_int inBondNum);

        template<typename _uni10_type>
          friend UniTensor<_uni10_type> permute( const UniTensor<_uni10_type>& T, uni10_int rowBondNum);

        template<typename _uni10_type>
          friend UniTensor<_uni10_type> partialTrace( const UniTensor<_uni10_type>& Tin, uni10_int la, uni10_int lb );

        template<typename _uni10_type>
          friend _uni10_type trace( const UniTensor<_uni10_type>& Tin );

        template<typename _uni10_type>
          friend void hosvd(UniTensor<_uni10_type>& Tin, uni10_int* group_labels, uni10_int* groups, uni10_uint64 groupsSize, 
              std::vector<UniTensor<_uni10_type> >& Us, UniTensor<_uni10_type>& S, std::vector<std::map<Qnum, Matrix<_uni10_type> > >& Ls, UNI10_INPLACE on);

        template<typename _uni10_type>
          friend void hosvd(UniTensor<_uni10_type>& Tin, uni10_int* group_labels, uni10_int* groups, uni10_uint64 groupsSize, 
              std::vector<UniTensor<_uni10_type> >& Us, UniTensor<_uni10_type>& S, std::vector<Matrix<_uni10_type> >& Ls, UNI10_INPLACE on);

        //bool similar(const UniTensor& Tb)const;
        //bool elemCmp(const UniTensor& UniT)const;
       
        template <typename _uni10_type>
          friend class UniTensor;

        template <typename _uni10_type>
          friend class Node;

        template <typename _uni10_type>
          friend class Network;

        friend class Network_dev;
        

    };

  template <typename uni10_type>
    uni10_int UniTensor<uni10_type>::COUNTER = 0;

  template <typename uni10_type>
    uni10_uint64 UniTensor<uni10_type>::ELEMNUM = 0;

  template <typename uni10_type>
    uni10_uint64 UniTensor<uni10_type>::MAXELEMNUM = 0;

  template <typename uni10_type>
    uni10_uint64 UniTensor<uni10_type>::MAXELEMTEN = 0;

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, const UniTensor<uni10_type>& UniT){

      if(UniT.paras == NULL){
          std::cout<<"This UniTensor has not been initialized."<< std::endl;
          return os;
      }

      if(!(*(UniT.status) & UniT.HAVEBOND)){
        if(UniT.U_elem->__ongpu)
          std::cout<<"\nScalar: " << UniT.U_elem->__elem[0]<<", onGPU";
        else{
          std::cout.precision(10);
          std::cout.setf(std::ios::fixed, std::ios::floatfield);
          std::cout<<"\nScalar: " << UniT.U_elem->__elem[0];
        }
        std::cout<<"\n\n";
        return os;
      }

      uni10_uint64 row = 0;
      uni10_uint64 col = 0;

      std::vector<Bond>bonds = (*UniT.bonds);
      for(uni10_uint64 i = 0; i < bonds.size(); i++)
        if(bonds[i].type() == BD_IN)
          row++;
        else
          col++;
      uni10_uint64 layer = std::max(row, col);
      uni10_uint64 nmlen = UniT.name->length() + 2;
      uni10_int star = 12 + (14 - nmlen) / 2;
      std::cout << std::endl;
      for(uni10_int s = 0; s < star; s++)
        std::cout << "*";
      if(UniT.name->length() > 0)
        std::cout << " " << (*UniT.name) << " ";
      for(uni10_int s = 0; s < star; s++)
        std::cout <<"*";
      std::cout<<std::endl;

      if(UniT.U_elem->__uni10_typeid == 1)
        std::cout << "REAL" << std::endl;
      else if(UniT.U_elem->__uni10_typeid == 2)
        std::cout << "COMPLEX" << std::endl;

      if(UniT.U_elem->__ongpu)
        std::cout<<"\n                 onGPU";
      std::cout << "\n             ____________\n";
      std::cout << "            |            |\n";
      uni10_uint64 llab = 0;
      uni10_uint64 rlab = 0;
      char buf[128];
      for(uni10_uint64 l = 0; l < layer; l++){
        if(l < row && l < col){
          llab = (*UniT.labels)[l];
          rlab = (*UniT.labels)[row + l];
          sprintf(buf, "    %5ld___|%-4d    %4d|___%-5ld\n", llab, bonds[l].dim(), bonds[row + l].dim(), rlab);
          std::cout<<buf;
        }
        else if(l < row){
          llab = (*UniT.labels)[l];
          sprintf(buf, "    %5ld___|%-4d    %4s|\n", llab, bonds[l].dim(), "");
          std::cout<<buf;
        }
        else if(l < col){
          rlab = (*UniT.labels)[row + l];
          sprintf(buf, "    %5s   |%4s    %4d|___%-5ld\n", "", "", bonds[row + l].dim(), rlab);
          std::cout << buf;
        }
        std::cout << "            |            |   \n";
      }
      std::cout << "            |____________|\n";

      std::cout << "\n================BONDS===============\n";
      for(uni10_uint64 b = 0; b < bonds.size(); b++)
        std::cout << bonds[b];

      std::cout<<"\n===============BLOCKS===============\n";
      typename std::map<Qnum, Block<uni10_type> >::const_iterator it = UniT.blocks->begin();
      for (; it != UniT.blocks->end(); it++ ){
        std::cout << "--- " << it->first << ": ";// << Rnum << " x " << Cnum << " = " << Rnum * Cnum << " ---\n\n";
        if((*(UniT.status) & UniT.HAVEELEM))
          std::cout<<it->second;
        else
          std::cout<<it->second.row() << " x "<<it->second.col()<<": "<<it->second.elemNum()<<std::endl<<std::endl;
      }
      std::cout << "Total elemNum: "<<(*UniT.U_elemNum)<<std::endl;
      std::cout << "====================================" << std::endl;
      os << "\n";
      return os;
    }

  template<typename uni10_type>
    UniTensor<uni10_type> operator*(const UniTensor<uni10_type>& Ta, uni10_double64 a){
      uni10_error_msg(!((*Ta.status) & Ta.HAVEELEM), "%s", "Cannot perform scalar multiplication on a tensor before setting its elements.");
      UniTensor<uni10_type> Tb(Ta);
      vectorScal(&a, Tb.U_elem, &Tb.U_elem->__elemNum);
      return Tb;
    }

  template<typename uni10_type>
    UniTensor<uni10_complex128> operator*(const UniTensor<uni10_type>& Ta, uni10_complex128 a){
      uni10_error_msg(!((*Ta.status) & Ta.HAVEELEM), "%s", "Cannot perform scalar multiplication on a tensor before setting its elements.");
      UniTensor<uni10_complex128> Tb(Ta);
      vectorScal(&a, Tb.U_elem, &Tb.U_elem->__elemNum);
      return Tb;
    }

  template<typename uni10_type>
    UniTensor<uni10_type> operator*(uni10_double64 a, const UniTensor<uni10_type>& Ta){return Ta*a;}

  template<typename uni10_type>
    UniTensor<uni10_complex128> operator*(uni10_complex128 a, const UniTensor<uni10_type>& Ta){return Ta*a;}

  template<typename T>
    UniTensor<T> operator+(const UniTensor<uni10_double64>& t1, const UniTensor<T>& t2){
      
      uni10_error_msg(t1.bondNum() != t2.bondNum(), "%s", "Cannot perform addition of two tensors having different bonds" );
      uni10_error_msg(!(*t1.bonds == *t2.bonds), "%s", "Cannot perform addition of two tensors having different bonds.");

      UniTensor<T> t3(t1);
      vectorAdd(t3.U_elem, t2.U_elem, &t3.U_elem->__elemNum);
      return t3;

    };

  template<typename T>
    UniTensor<uni10_complex128> operator+(const UniTensor<uni10_complex128>& t1, const UniTensor<T>& t2){
      
      uni10_error_msg(t1.bondNum() != t2.bondNum(), "%s", "Cannot perform addition of two tensors having different bonds" );
      uni10_error_msg(!(*t1.bonds == *t2.bonds), "%s", "Cannot perform addition of two tensors having different bonds.");

      UniTensor<uni10_complex128> t3(t1);
      vectorAdd(t3.U_elem, t2.U_elem, &t3.U_elem->__elemNum);
      return t3;

    };

  template<typename T>
    UniTensor<T> operator-(const UniTensor<uni10_double64>& t1, const UniTensor<T>& t2){
      
      uni10_error_msg(t1.bondNum() != t2.bondNum(), "%s", "Cannot perform addition of two tensors having different bonds" );
      uni10_error_msg(!(*t1.bonds == *t2.bonds), "%s", "Cannot perform addition of two tensors having different bonds.");

      UniTensor<T> t3(t1);
      vectorSub(t3.U_elem, t2.U_elem, &t3.U_elem->__elemNum);
      return t3;

    };

  template<typename T>
    UniTensor<uni10_complex128> operator-(const UniTensor<uni10_complex128>& t1, const UniTensor<T>& t2){
      
      uni10_error_msg(t1.bondNum() != t2.bondNum(), "%s", "Cannot perform addition of two tensors having different bonds" );
      uni10_error_msg(!(*t1.bonds == *t2.bonds), "%s", "Cannot perform addition of two tensors having different bonds.");

      UniTensor<uni10_complex128> t3(t1);
      vectorSub(t3.U_elem, t2.U_elem, &t3.U_elem->__elemNum);
      return t3;

    };

};

#endif
