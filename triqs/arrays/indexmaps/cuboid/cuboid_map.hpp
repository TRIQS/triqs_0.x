
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef TRIQS_ARRAYS_INDEXMAP_CUBOIDMAP_H
#define TRIQS_ARRAYS_INDEXMAP_CUBOIDMAP_H

#include "./cuboid_domain.hpp"
#include "./index_order.hpp"
#include "./cuboid_map_iterator.hpp"
#include <vector>
#include <boost/preprocessor/repetition.hpp>

namespace triqs { namespace arrays {

 namespace Tag { struct cuboid_map:indexmap {}; };
 namespace indexmaps { 

  template< bool BC, class D, class K> struct _chk;
  template< class D, class K> struct _chk<true, D, K>  { static void invoke (D const & d, K const & k) {d.assert_key_in_domain(k);} };
  template< class D, class K> struct _chk<false, D, K> { static void invoke (D const & d, K const & k) {} };

  /** Standard hyper_rectangular arrays, implementing the IndexMap concept.  */
  template<typename IndexOrderType, bool CheckBounds>
   class cuboid_map : public Tag::cuboid_map { 
    public : 

     typedef cuboid_domain< IndexOrderType::rank> domain_type;
     domain_type const & domain() const { return mydomain;}

     cuboid_map ():mydomain(), start_shift_(0){}
     cuboid_map(domain_type const & C):mydomain(C), start_shift_(0){compute_stride_compact();} 
     cuboid_map (cuboid_map const & C):mydomain(C.mydomain), strides_(C.strides_), start_shift_(C.start_shift_) {}

     /// Returns the shift in position of the element key. 
     template <typename KeyType>
      size_t operator[] (KeyType const & key ) const {
       _chk<CheckBounds, domain_type, KeyType>::invoke (this->domain(),key);
       return start_shift_ + dot_product(key,this->strides());
      }

     typedef cuboid_map_iterator<cuboid_map, typename IndexOrder::OptimalTraversal<IndexOrderType>::type > iterator; 

     friend std::ostream & operator << (std::ostream & out, const cuboid_map & P) {
      return out <<"  ordering = {" << typename index_order_type::perm()<<"}"<<std::endl
       <<"  Lengths  :  "<<P.lengths() << std::endl 
       <<"  Stride  : "<<P.strides_ << std::endl; 
     }

     //------------- end of concept  ---------------------

     static const unsigned int rank = IndexOrderType::rank;
     typedef mini_vector<size_t,rank> lengths_type;
     typedef mini_vector<std::ptrdiff_t, rank> strides_type;
     typedef IndexOrderType index_order_type;

     /// Construction from the length, the stride, start_shift 
     cuboid_map(lengths_type const & Lengths, strides_type const & strides, std::ptrdiff_t start_shift ): 
      mydomain(Lengths), strides_(strides), start_shift_(start_shift){ }

     /// Construction from another cuboid_map with the same order (used in grouping indices)
     template<typename IO> cuboid_map (cuboid_map<IO,CheckBounds> const & C):
      mydomain(C.domain()), strides_(C.strides()), start_shift_(C.start_shift()) {
      static_assert( (Permutations::is_equal<typename IO::perm,typename IndexOrderType::perm>::value), "Internal error");
      }

     /// 
     bool is_contiguous() const {   
      const size_t last_index = IndexOrderType::template memory_rank_to_index<rank-1>::value; 
      return (strides()[last_index] * this->lengths()[last_index] == mydomain.number_of_elements()); 
     }

     /// An iterator in a "natural" order ... 
     typedef cuboid_map_iterator<cuboid_map, Permutations::identity<rank> > natural_iterator; 

     IndexOrderType index_order () const {return IndexOrderType();}

     size_t start_shift() const { return start_shift_;}

     lengths_type const & lengths() const { return mydomain.lengths();}

     strides_type const & strides() const { return this->strides_;}

    protected:
     domain_type mydomain;
     strides_type strides_;
     std::ptrdiff_t start_shift_;

    private :
     //  BOOST Serialization
     friend class boost::serialization::access;
     template<class Archive> void serialize(Archive & ar, const unsigned int version) {
      ar & boost::serialization::make_nvp("domain",mydomain);
      ar & boost::serialization::make_nvp("strides",strides_);
      ar & boost::serialization::make_nvp("start_shift",start_shift_);
     }
     void compute_stride_compact();// impl below

   }; //------------- end class ---------------------

  template<typename A1, typename A2, bool B1, bool B2>
   bool compatible_for_assignment (const cuboid_map<A1,B1> & X1, const cuboid_map<A2,B2> & X2) { return X1.lengths() == X2.lengths();}

  template<typename A1, typename A2, bool B1, bool B2>
   bool raw_copy_possible (const cuboid_map<A1,B1> & X1,const cuboid_map<A2,B2> & X2) {
    return ( (X1.index_order() == X2.index_order())  
      && X1.is_contiguous() && X2.is_contiguous() 
      && (X1.domain().number_of_elements()==X2.domain().number_of_elements()));
   }

  template<typename P, typename I, typename S2, bool B2>
   struct indexmap_iterator_adapter< cuboid_map_iterator<P,I> , cuboid_map<S2,B2> > {
    typedef cuboid_map_iterator<cuboid_map<S2,B2> , typename cuboid_map_iterator<P,I>::iteration_order>  type;
   };

  /*-----------------------------------------
   *       Implementation details 
   *----------------------------------------*/

  namespace cuboid_details {

   // go over the indices in the order of their storage in memory
   template<typename IndexOrderType, typename T1, typename T2, int rank, int v> struct csc_impl {
    static const size_t u =  IndexOrderType::template memory_rank_to_index <rank-v>::value; 
    static void do_it( T1 const & lengths_, T2 & strides_, size_t & str) {
     strides_[u]  = str; 
     str *= lengths_ [u]; 
     csc_impl<IndexOrderType,T1,T2,rank, v-1>::do_it(lengths_,strides_,str);
    }
   };

   template<typename IndexOrderType, typename T1, typename T2, int rank> 
    struct csc_impl<IndexOrderType,T1,T2,rank,0> { static void do_it( T1 const & , T2 & , size_t & ) {} };
  }

  template<typename IndexOrderType, bool BC>
   void cuboid_map<IndexOrderType, BC>::compute_stride_compact() {
    size_t str = 1;
    cuboid_details::csc_impl<index_order_type,lengths_type, strides_type,rank,rank>::do_it(this->lengths(),this->strides_, str);
    assert(this->domain().number_of_elements()==str);
   }

 }
}}//namespace triqs::arrays 
#endif
