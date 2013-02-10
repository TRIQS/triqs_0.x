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
#ifndef TRIQS_ARRAYS_INDEXMAP_CUBOID_MAP_H
#define TRIQS_ARRAYS_INDEXMAP_CUBOID_MAP_H
#include "./domain.hpp"
#include "./mem_layout.hpp"
#include "../../impl/flags.hpp"
#include <vector>
#include <boost/iterator/iterator_facade.hpp>

namespace triqs { namespace arrays { namespace indexmaps { namespace cuboid { 

 template< bool BC, class D, class K> struct _chk;
 template< class D, class K> struct _chk<true, D, K>  { static void invoke (D const & d, K const & k) {d.assert_key_in_domain(k);} };
 template< class D, class K> struct _chk<false, D, K> { static void invoke (D const & d, K const & k) {} };

 template< bool BC, class D, class ... K> struct _chk_v;
 template< class D, class ... K> struct _chk_v<true, D, K...>  { static void invoke (D const & d, K const & ... k) {d.assert_key_in_domain_v(k...);} };
 template< class D, class ... K> struct _chk_v<false, D, K...> { static void invoke (D const & d, K const & ... k) {} };


 /** Standard hyper_rectangular arrays, implementing the IndexMap concept.  
 */
 template<int Rank, ull_t OptionsFlags>
  class map {
   static constexpr bool CheckBounds = flags::bound_check(OptionsFlags);
   public : 

   typedef domain<Rank> domain_type;
   domain_type const & domain() const { return mydomain;}

   //map (memory_layout<Rank> ml):mydomain(), start_shift_(0), memory_order_(ml) {} 
   map (memory_layout<Rank> ml = memory_layout<Rank>()):mydomain(), start_shift_(0), memory_order_(ml) {} 

   map (map const & C) = default; //:mydomain(C.mydomain), strides_(C.strides_), start_shift_(C.start_shift_), memory_order_  {}

   map(domain_type const & C):
    mydomain(C), start_shift_(0), memory_order_( mem_layout::c_order(Rank) ) { compute_stride_compact(); } 

   map(domain_type const & C, memory_layout<Rank> ml):
    mydomain(C), start_shift_(0), memory_order_(ml) { compute_stride_compact(); } 

   /// Returns the shift in position of the element key. 
   template <typename KeyType>
    size_t operator[] (KeyType const & key ) const {
     _chk<CheckBounds, domain_type, KeyType>::invoke (this->domain(),key);
     return start_shift_ + dot_product(key,this->strides());
    }

   friend std::ostream & operator << (std::ostream & out, const map & P) {
    out <<"  ordering = {";
    //arrays::permutations::print(out, P.memory_order_); 
    return out <<"}"<<std::endl
     <<"  Lengths  :  "<<P.lengths() << std::endl 
     <<"  Stride  : "<<P.strides_ << std::endl; 
   }

   //------------- end of concept  ---------------------

   // mv into the concept ?
   template<typename ... Args> size_t operator()(Args const & ... args) const { 
     _chk_v<CheckBounds, domain_type, Args...>::invoke (this->domain(),args...);
    return start_shift_ + _call_impl<0>(args...);
   }
   private : 
   template<int N, typename Arg0, typename ... Args> 
    size_t _call_impl( Arg0 const & arg0, Args const & ... args) const { return arg0* strides_[N] + _call_impl<N+1>(args...); }
   template<int N> size_t _call_impl() const { return 0;}

   public:
   static const unsigned int rank = Rank; //IndexOrderType::rank;
   typedef mini_vector<size_t,rank> lengths_type;
   typedef mini_vector<std::ptrdiff_t, rank> strides_type;

   /// Construction from the length, the stride, start_shift 
   map(lengths_type const & Lengths, strides_type const & strides, std::ptrdiff_t start_shift ): 
    mydomain(Lengths), strides_(strides), start_shift_(start_shift){ }

   /*/// Construction from another map with the same order (used in grouping indices)
     template<typename IO> map (map<IO,CheckBounds> const & C):
     mydomain(C.domain()), strides_(C.strides()), start_shift_(C.start_shift()) {
     static_assert( (Permutations::is_equal<typename IO::perm,typename IndexOrderType::perm>::value), "Internal error");
     }
     */
   /// 
   bool is_contiguous() const {   
    const size_t last_index = mem_layout::memory_rank_to_index(memory_order_.value, rank-1);
    //const size_t last_index = IndexOrderType::template memory_rank_to_index<rank-1>::value; 
    return (strides()[last_index] * this->lengths()[last_index] == mydomain.number_of_elements()); 
   }

   //IndexOrderType index_order () const {return IndexOrderType();}
   size_t start_shift() const { return start_shift_;}
   lengths_type const & lengths() const { return mydomain.lengths();}
   strides_type const & strides() const { return this->strides_;}

   memory_layout<Rank> memory_indices_layout() const { return memory_order_;}
    
   bool memory_layout_is_c() const { return memory_indices_layout() == memory_layout<Rank> (mem_layout::c_order(Rank));}
   bool memory_layout_is_fortran() const { return memory_indices_layout() == memory_layout<Rank> (mem_layout::fortran_order(Rank));}

   private :
   domain_type mydomain;
   strides_type strides_;
   std::ptrdiff_t start_shift_;
   memory_layout<Rank> memory_order_;

   //  BOOST Serialization
   friend class boost::serialization::access;
   template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::make_nvp("domain",mydomain);
    ar & boost::serialization::make_nvp("strides",strides_);
    ar & boost::serialization::make_nvp("start_shift",start_shift_);
   }
   // for construction
   void compute_stride_compact() { 
    size_t str = 1;
    csc_impl(memory_order_.value, str, std::integral_constant<int,rank>());
    assert(this->domain().number_of_elements()==str);
   }
   template<int v>
    void csc_impl(ull_t order, size_t & str, std::integral_constant<int,v>) {
     size_t u = mem_layout::memory_rank_to_index(order, rank-v);
     this->strides_[u]  = str; 
     str *= this->lengths() [u]; 
     csc_impl(order, str, std::integral_constant<int,v-1>());
    }
   void csc_impl(ull_t order,size_t & str, std::integral_constant<int,0>) {}

   public:

   /**
    * Iterator on a cuboid_map, modeling the IndexMapIterator concept.
    * Iteration order is the order in which to iterate on the indices.
    * It is given by a permutation, with the same convention as IndexOrder.
    */
   class iterator : public boost::iterator_facade< iterator, const std::ptrdiff_t, boost::forward_traversal_tag > {
    public:
     typedef map indexmap_type;
     typedef typename domain_type::index_value_type indices_type;
     typedef const std::ptrdiff_t return_type;
     iterator (): io(0), im(NULL), pos(0),atend(true) {}
     iterator (const map & P, bool atEnd=false, ull_t iteration_order=0): 
      io((iteration_order !=0 ? iteration_order : P.memory_indices_layout().value)),
      im(&P), pos(im->start_shift()),atend(atEnd) {}
     indices_type const & indices() const { return indices_tuple; }
     operator bool() const { return !atend;}
    private:
     friend class boost::iterator_core_access;
     void increment(){ inc_ind_impl (std::integral_constant<int,Rank>()); }
     template<int v> void inc_ind_impl(std::integral_constant<int,v>) { 
      const size_t p = permutations::apply(io, rank-v);
#ifdef TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
      if (atend) TRIQS_RUNTIME_ERROR << "Iterator in cuboid can not be pushed after end !";
#endif
      if (indices_tuple[p] < im->lengths()[p]-1) { ++(indices_tuple[p]); pos += im->strides()[p]; return; }
      indices_tuple[p] = 0; 
      pos -= (im->lengths()[p]-1) * im->strides()[p];
      inc_ind_impl (std::integral_constant<int,v-1>());
     }
     void inc_ind_impl(std::integral_constant<int,0>) { atend = true;} 
     bool equal(iterator const & other) const {return ((other.im==im)&&(other.atend==atend)&&(other.pos==pos));}
     return_type & dereference() const { assert (!atend); return pos; }
     ull_t io; 
     map const * im;
     indices_type indices_tuple; 
     std::ptrdiff_t pos;
     bool atend;
   };

  }; //------------- end class ---------------------
}//namespace cuboid

template<int R1, int R2, ull_t OptFlags1, ull_t OptFlags2>
bool compatible_for_assignment (const cuboid::map<R1,OptFlags1> & X1, const cuboid::map<R2,OptFlags2> & X2) { return X1.lengths() == X2.lengths();}

template<int R1, int R2, ull_t OptFlags1, ull_t OptFlags2>
 bool raw_copy_possible (const cuboid::map<R1,OptFlags1> & X1,const cuboid::map<R2,OptFlags2> & X2) {
  return ( (X1.memory_indices_layout() == X2.memory_indices_layout())  
    && X1.is_contiguous() && X2.is_contiguous() 
    && (X1.domain().number_of_elements()==X2.domain().number_of_elements()));
 }

/*
/// TO BE REPLACED AT RUN TIME 
template<typename I, int R2, bool B2>
struct indexmap_iterator_adapter< I , cuboid::map<R2,B2> > { typedef typename cuboid::map<R2,B2>::iterator type; };
*/
}}}//namespace triqs::arrays::indexmaps
#endif
