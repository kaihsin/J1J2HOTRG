/****************************************************************************
 *  @file Network.h
 *  @license
 *    Universal Tensor Network Library
 *    Copyright (c) 2013-2014
 *    National Taiwan University
 *    National Tsing-Hua University

 *
 *    This file is part of Uni10, the Universal Tensor Network Library.
 *
 *    Uni10 is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU Lesser General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    Uni10 is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public License
 *    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
 *  @endlicense
 *  @brief Header file for Newtork class
 *  @author Yun-Da Hsieh
 *  @date 2014-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/
#ifndef NETWORK_H
#define NETWORK_H

#include "uni10/uni10_api/UniTensor.h"
#include "uni10/uni10_api/network_tools/Node.h"
#include "uni10/uni10_api/network_tools/network_tools.h"

namespace uni10 {

  template<typename uni10_type>
    class Network;

  template<typename uni10_type>
    std::ostream& operator<< (std::ostream& os, Network<uni10_type>& _b);

  ///@class Network
  ///@brief The Network class defines the tensor networks
  ///
  /// A Network consists of connections which is specified by labels.
  /// To construct a network, prepare a file as the following example,
  ///
  ///     A: 1 2; 3 4
  ///     B: 3 4; 5 6
  ///     C: -1; 7 1
  ///     D: -2; 2 8
  ///     E: 7 5; -3
  ///     F: 6 8 -4
  ///     TOUT: -1 -2; -3 -4
  ///     ORDER: A B C E D F
  ///
  /// The first column specifies the labels of the tensors. The line starting with `TOUT` specifies the labels
  /// of the output tensor.
  /// Labels are separated by space or comma.  Incoming and outgoing labels are seperated by a semi-colon.
  /// The last line, starting from “ORDER: ” is optional. This line is used to provide a suggested order to
  /// construct the pair-wise contractions. The contraction order can be forced by putting parentheses at
  /// appropriated place as,
  ///
  ///     ORDER: ((((A B) C) E) (D F))
  ///
  /// @note The `TOUT:` line is required. If the result is a scalar, keep the line `TOUT:` without any labels.
  ///
  /// @see UniTensor
  /// @example egN1.cpp

  template <typename uni10_type>
    class Network {
      public:
        /// @brief Construct Network
        ///
        /// Constructs a Network, initializing with the file \c fname. The file specifies the connections between
        /// tensors in the network.
        /// @param fname %Network filename
        Network(const std::string& fname);

        /// @brief Construct Network
        ///
        /// Constructs with the file fname and put the UniTensor's \c tens into the network.
        /// The UniTensor's in the array \c tens are assigned in the order as in the network file \c fname.
        /// @param fname %Network filename
        /// @param tens Array of tensors
        Network(const std::string& fname, const std::vector<UniTensor<uni10_type>*>& tens);

        /// @brief Destructor
        ///
        /// Destroys Network and frees all the intermediate tensors.
        ~Network();

        /// @brief Assign tensor to Network
        ///
        /// Assigns UniTensor \c uT to position \c idx in Network
        /// If \c force is set \c true, replace the tensor without reconstructing the pair-wise contraction sequence.
        /// @param idx Position
        /// @param UniT A UniTensor
        /// @param force If set \true, replace without chaning the contraction sequence. Defaults to \c true.
        void putTensor(size_t idx, const UniTensor<uni10_type>& UniT, bool force=true);
        /// @overload
        void putTensor(size_t idx, const UniTensor<uni10_type>* UniT, bool force=true);

        /// @brief Assign tensor to Network
        ///
        /// Assigns UniTensor \c uT to the position labeled by \c name in Network
        /// If \c force is set \c true, replace the tensor without reconstructing the pair-wise contraction sequence.
        /// @param name Name of tensor in Network
        /// @param UniT A UniTensor
        /// @param force If set \true, replace without chaning the contraction sequence. Defaults to \c true.
        void putTensor(const std::string& name, const UniTensor<uni10_type>& UniT, bool force=true);
        /// @overload
        void putTensor(const std::string& name, const UniTensor<uni10_type>* UniT, bool force=true);

        /// @brief Assign tensor to  Network
        ///
        /// Assigns the transpose of \c uT  to the position labeled by \c nameT in Network.
        /// If \c force is set \c true, replace the tensor without reconstructing the pair-wise contraction sequence.
        /// @param nameT Name of tensor in Network
        /// @param UniT A UniTensor
        /// @param force If set \true, replace without chaning the contraction sequence. Defaults to \c true.
        void putTensorT(const std::string& nameT, const UniTensor<uni10_type>& UniT, bool force=true);
        /// @overload
        void putTensorT(const std::string& nameT, const UniTensor<uni10_type>* UniT, bool force=true);

        /// @brief Contract Network
        ///
        /// Performs contraction of tensors in Network, returns a UniTensor named \c name.
        /// @param name Name of the result tensor
        /// @return A UniTensor
        UniTensor<uni10_type> launch(const std::string& name="");

        void launch(UniTensor<uni10_type>& Tout, const std::string& name="");
        /// @brief Print out the memory usage
        /// Prints out the memory usage and requirement to contract  Network as:
        /** @code
          ===== Network profile =====
          Memory Requirement: 1032
          Sum of memory usage: 504
          Maximun tensor:
          elemNum: 19
          4 bonds and labels: 1, 2, 3, 4,
          ===========================
          @endcode
         **/
        /// In the above example, to contract Network, the memory requirement is 1032 bytes.
        /// The maximum tensor in Network has 19 elements and has four bonds with labels 1, 2, 3, 4.
        std::string profile(bool print=true);
        /// @brief Print out Network
        ///
        /// For a newtork described in the following network file,
        /** @code
W1: -1; 0 1 3
W2: -2; 7 10 11
U: 3 7; 4 8
Ob: 1 4; 2 5
UT: 5 8; 6 9
W1T: 0 2 6; -3
W2T: 9 10 11; -4
Rho: -3 -4; -1 -2
TOUT:
ORDER: W1 W1T W2 W2T U Ob UT Rho
@endcode */
        /// Before calling \ref launch(), the function prints out  Network as
        ///
        /** @code
W1: i[-1] o[0, 1, 3]
W2: i[-2] o[7, 10, 11]
U: i[3, 7] o[4, 8]
Ob: i[1, 4] o[2, 5]
UT: i[5, 8] o[6, 9]
W1T: i[0, 2, 6] o[-3]
W2T: i[9, 10, 11] o[-4]
Rho: i[-3, -4] o[-1, -2]
TOUT:
@endcode
         */
        ///
        /// `i` are labels for incoming bonds and `o` for outgoing bonds.
        ///  After launch() is called, in addition to the network, it also prints out the contraction sequence, or
        /// the binary tree of the pair-wise contractions is also printed as,
        ///
        /** @code
         *(1):
         |   *(70): 7, 9, -4, -2,
         |   |   *(70): -1, 7, 9, -3,
         |   |   |   *(70): -1, 0, 7, 2, 6, 9,
         |   |   |   |   W1(20): -1, 0, 1, 3,
         |   |   |   |   *(20): 3, 7, 1, 2, 6, 9,
         |   |   |   |   |   *(20): 3, 7, 8, 1, 2, 5,
         |   |   |   |   |   |   U(6): 3, 7, 4, 8,
         |   |   |   |   |   |   Ob(6): 1, 4, 2, 5,
         |   |   |   |   |   UT(6): 5, 8, 6, 9,
         |   |   |   W1T(20): 0, 2, 6, -3,
         |   |   Rho(924): -3, -4, -1, -2,
         |   *(70): -2, 7, 9, -4,
         |   |   W2(20): -2, 7, 10, 11,
         |   |   W2T(20): 9, 10, 11, -4,
         @endcode
         */
        ///
        /// The output shows how the network is contracted.

        friend std::ostream& operator<< <>(std::ostream& os, Network& net);

      private:

        std::string fname;
        void preprint(std::ostream& os, Node<uni10_type>* nd, uni10_int32 layer)const;  //pre-order print
        std::vector<std::string> names;
        std::map<std::string, uni10_uint64> name2pos;
        std::vector< std::vector<uni10_int32> > label_arr;
        std::vector< uni10_int32 > Rnums;
        std::vector<Node<uni10_type>*> leafs;
        std::vector<UniTensor<uni10_type>*> tensors;
        std::vector< std::vector<_Swap> > swaps_arr;
        std::vector<uni10_bool> swapflags;
        std::vector<uni10_int32> conOrder;  //contraction order;
        std::vector<uni10_int32> order; //add order
        std::vector<uni10_int32> brakets;   //add order
        Node<uni10_type>* root;
        bool load;  //whether or not the network is ready for contraction, construct=> load=true, destruct=>load=false
        int times;  //construction times
        int tot_elem;   //total memory ussage
        int max_elem;   //maximum
        void construct();
        void destruct();
        void matching(Node<uni10_type>* sbj, Node<uni10_type>* tar);
        void branch(Node<uni10_type>* sbj, Node<uni10_type>* tar);
        UniTensor<uni10_type> merge(Node<uni10_type>* nd);
        void clean(Node<uni10_type>* nd);
        void fromfile(const std::string& fname);
        void findConOrd(Node<uni10_type>* nd);
        void addSwap();
        int rollcall();
        size_t sum_of_memory_usage();
        size_t max_tensor_elemNum();
        size_t memory_requirement();
        void _max_tensor_elemNum(Node<uni10_type>* nd, uni10_uint64& max_num, Node<uni10_type>& max_nd) const;
        size_t _sum_of_tensor_elem(Node<uni10_type>* nd) const;
        size_t _elem_usage(Node<uni10_type>* nd, uni10_uint64& usage, uni10_uint64& max_usage)const;
        bool hasOrder;
        int order_pos;
    };

  template <typename uni10_type>
    std::ostream& operator<< (std::ostream& os, Network<uni10_type>& net){
      os<<std::endl;
      for(int i = 0; i < net.names.size(); i++){
        os<<net.names[i]<< ": ";
        if(net.Rnums[i])
          os<<"i[";
        for(int l = 0; l < net.Rnums[i]; l++){
          os<<net.label_arr[i][l];
          if(l < net.Rnums[i] - 1)
            os<<", ";
        }
        if(net.Rnums[i])
          os<<"] ";
        if(net.label_arr[i].size() - net.Rnums[i])
          os<<"o[";
        for(int l = net.Rnums[i]; l < net.label_arr[i].size(); l++){
          os<<net.label_arr[i][l];
          if(l < net.label_arr[i].size() - 1)
            os<<", ";
        }
        if(net.label_arr[i].size() - net.Rnums[i])
          os<<"]";
        os<<std::endl;
      }
      os<<std::endl;
      if(net.rollcall() < 0)
        net.preprint(os, net.root, 0);
      else
        os<<"\nSome tensors have not yet been given!\n\n";
      return os;
    }

};  /* namespace uni10 */
#endif /* NETWORK<uni10_type>_H */
