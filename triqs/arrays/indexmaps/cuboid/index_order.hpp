
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
#ifndef TRIQS_ARRAYS_INDEXMAP_STORAGE_ORDER_H
#define TRIQS_ARRAYS_INDEXMAP_STORAGE_ORDER_H

#include "../permutation.hpp"
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/vector.hpp>
#include "../../impl/tuple_tools.hpp"

namespace triqs { namespace arrays { 

 namespace Tag { struct index_order_static{};}
 namespace indexmaps { namespace IndexOrder { 

  /*
   * index_to_memory_rank  : i ---> r : to the index (0,1, ... etc), associates the rank in memory
   *                         e.g. r=0 : fastest index, r = RANK-1 : the slowest
   * memory_rank_to_index  : the inverse mapping : r---> i : 
   *                         0-> index of the fastest index, etc..
   */

  /// Storage using C conventions
  template <int RANK> 
   struct C: Tag::index_order_static {
    static const unsigned int rank = RANK;
    template<int r> struct index_to_memory_rank { static const int value =  (RANK - 1) - r;};
    template<int r> struct memory_rank_to_index { static const int value =  (RANK - 1) - r;};
    typedef Permutations::reverse_identity<rank> perm; 
    static bool is_Fortran() { return (rank <=1);}
    static bool is_C() { return true;}
    template<class Archive> void serialize(Archive & ar, const unsigned int version) {}
    static const char C_or_F = 'C';
   };

  /// Storage using Fortran conventions
  template <int RANK> 
   struct Fortran: Tag::index_order_static  {
    static const unsigned int rank = RANK;
    template<int r> struct index_to_memory_rank { static const int value =  r;};
    template<int r> struct memory_rank_to_index { static const int value =  r;};
    typedef Permutations::identity<rank> perm; 
    static bool is_Fortran()  { return true;}
    static bool is_C()  { return (rank<=1);}
    template<class Archive> void serialize(Archive & ar, const unsigned int version) {}
    static const char C_or_F = 'F';
   };

  /** 
   * Storage using custom order conventions
   * - P must be a Permutation. 
   *   if P[0] : the fastest index, 
   *   if P[RANK-1] : the slowest index
   *   Example : 
   *   custom<Permutation<0,1,2> >  : same as index_order_Fortran, the first index is the fastest
   *   custom<Permutation<2,1,0> >  : same as index_order_C, the last index is the fastest
   *   custom<Permutation<1,2,0> >  : storage (i,j,k) is : index j is fastest, then k, then i 
   *     
   */
  template <typename P> 
   struct custom: Tag::index_order_static  {
    typedef P perm_inv; 
    typedef Permutations::inverse<P> perm;
    static const unsigned int rank = Permutations::length<P>::value;
    template<int r> struct index_to_memory_rank { static const int value =  Permutations::eval<perm, r >::value; };
    template<int r> struct memory_rank_to_index { static const int value =  Permutations::eval<perm_inv, r >::value; };
    static bool is_Fortran() { return Permutations::is_equal<P,typename Fortran<rank>::perm >::value; }
    static bool is_C() { return Permutations::is_equal<P,typename C<rank>::perm >::value; }
    template<class Archive> void serialize(Archive & ar, const unsigned int version) {}
    static const char C_or_F = 'N';
   };

  //-----------------------------------------------------------

  template<typename ORDER> struct OptimalTraversal {};

  template<int RANK> struct OptimalTraversal<C<RANK> > { typedef Permutations::reverse_identity<RANK> type; };
  template<int RANK> struct OptimalTraversal<Fortran<RANK> > { typedef Permutations::identity<RANK> type; };
  template<typename P> struct OptimalTraversal<custom<P> > { typedef P type; };

  //-----------------------------------------------------------

  template<typename MemoryOrderTYPE,typename ArgsTuple> struct SlicedIndexOrder;

  template<int RANK, typename ArgsTuple>
   struct SlicedIndexOrder<C<RANK>, ArgsTuple> {
    static const int NN=RANK - TupleTools::CountHowManyInt<ArgsTuple>::value;
    typedef C<NN> type;
   };

  template<int RANK,typename ArgsTuple>
   struct SlicedIndexOrder<Fortran<RANK>, ArgsTuple> {
    static const int NN=RANK - TupleTools::CountHowManyInt<ArgsTuple>::value;
    // Changg Fortran<1> into C<1>. Fortran<1> is never returned
    typedef typename boost::mpl::if_c< NN==1, C<1>, Fortran<NN> >::type type;
   };

  /*
   * Given a general IndexOrder given by a permutation P, we want to compute the 
   * IndexOrder resulting from the slicing with arguments (A0,....,A_r-1)
   *
   * Done in two steps : 
   *  1) Compute a mask : (A0,...,A_r-1) is replaced by integers, where
   *      * if A0 is a range : -1
   *      * else : a increasing integer
   *    e.g. (4,Range(xxx), Range(xxx), 5) ----> (-1,0,1,-1).
   *
   *    Formula : 
   *      mask(n,p) = [p if A[n] is range else -1] + mask (n+1,p + (1 if A[n] is range else 0))
   *      mask(r,p ) = []
   *
   *  2) Compute the new permutation determining the sliced IndexOrder, from P, the mask M= mask(0,0)     
   *     p(n) = [ M[P[n]] ] + p(n+1)  if   M[P[n]]!=-1 else  p(n+1)
   *     p(r) = [] 
   *     
   */
  namespace details {

   template <int N, int P, typename ArgsTuple, int c> struct mask {
    static const bool is_range = boost::is_base_of<typename boost::tuples::element<N,ArgsTuple>::type, range >::type::value ;
    static const int head = boost::mpl::if_c<is_range, boost::mpl::int_<P>, boost::mpl::int_<-1> >::type::value;
    typedef Permutations::push_front<head, typename mask<N+1,P + is_range,ArgsTuple,c-1>::type > type;  
   };
   template <int N, int P, typename ArgsTuple> struct mask<N,P,ArgsTuple,0> { typedef  Permutations::empty type; };

   template<int N, typename ArgsTuple, typename perm_inv, typename Mask, int c > struct slo_impl {
    static const int pn = Permutations::eval<perm_inv,N>::value;
    static const int r2 = Permutations::eval<Mask,pn>::value;
    typedef typename slo_impl<N+1, ArgsTuple,perm_inv,Mask,c-1>::type previous;
    typedef typename boost::mpl::if_c<(r2 != -1), Permutations::push_front<r2,previous>, previous>::type type;
   };
   template<int N, typename A, typename pinv, typename M > struct slo_impl<N,A,pinv,M,0> { typedef Permutations::empty type; };

  }//details

  template<typename P,typename ArgsTuple>
   struct SlicedIndexOrder<custom<P>, ArgsTuple> { 
    typedef typename details::mask<0,0,ArgsTuple,custom<P>::rank >::type M;
    typedef typename  details::slo_impl<0,ArgsTuple,typename custom<P>::perm_inv, M, custom<P>::rank >::type P2;
    typedef custom<P2> type;
   };
 }}
 //-----------------------------------------------------------
 // checks if two static order are identical
 template<typename T1, typename T2>
  typename boost::enable_if<boost::mpl::and_< triqs::arrays::Tag::check<triqs::arrays::Tag::index_order_static,T1>, 
	   triqs::arrays::Tag::check<triqs::arrays::Tag::index_order_static,T2> >, bool >::type 
	    operator ==( T1 const & x1, T2 const & x2)  { 
	     return ( (x1.rank==x2.rank) && Permutations::is_equal<typename T1::perm,typename T2::perm>::value) ;
	    }
}}//namespace triqs::arrays 
#endif
