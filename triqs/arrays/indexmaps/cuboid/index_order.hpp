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
#include "../permutation2.hpp"
namespace triqs { namespace arrays { namespace indexmaps { namespace index_order {

  /*
   * The storage order is given by a permutation stored in a ull_t (unsigned long long) as in permutations::..
   *
   * If P is this permutation :  
   *   P[0] : the fastest index, 
   *   P[RANK-1] : the slowest index
   *   Example : 
   *   012 : Fortran, the first index is the fastest
   *   210:  C the last index is the fastest
   *   120 : storage (i,j,k) is : index j is fastest, then k, then i 
   * 
   * index_to_memory_rank  : i ---> r : to the index (0,1, ... etc), associates the rank in memory
   *                         e.g. r=0 : fastest index, r = RANK-1 : the slowest
   * memory_rank_to_index  : the inverse mapping : r---> i : 
   *                         0-> the fastest index, etc..
   *
   * All these computations can be done *at compile time* (constexpr)
   */
  
  constexpr int memory_rank_to_index(ull_t p, int r) { return permutations::apply(p, r);} 
  constexpr int index_to_memory_rank(ull_t p, int r) { return permutations::apply(permutations::inverse(p), r);} 
  constexpr int rank(ull_t p) { return permutations::size(p);}
  constexpr ull_t optimal_traversal(ull_t p) { return permutations::inverse(p);}

  constexpr bool is_fortran (ull_t p){ return p == permutations::identity(permutations::size(p));}
  constexpr bool is_c       (ull_t p){ return p == permutations::ridentity(permutations::size(p));}

  constexpr ull_t fortran_order (int n){ return permutations::identity(n);}
  constexpr ull_t c_order       (int n){ return permutations::ridentity(n);}

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

  /* template <ull_t MemoryOrder, int P, int c, typename... SliceArgs> struct mask;
     template <ull_t MemoryOrder, int P, int c>                        struct mask<MemoryOrder,P,c> { static constexpr ull_t value = 0;};
     template <ull_t MemoryOrder, int P, int c, typename A0, typename... SliceArgs> struct mask<MemoryOrder,c,A0,SliceArgs...> { 
     static constexpr bool is_range = std::is_base_of<range, A0>::value; 
     static constexpr ull_t value = ((is_range ? P+1 : 0) << (4*c)) + mask<MemoryOrder,P + is_range, c+1,SliceArgs...>::value;
     };


     template <typename... SliceArgs> struct mask;
     template <>                        struct mask<P,c> { static constexpr ull_t invoke(ull_t mo, int P, int c) { return 0ull;}};
     template <typename A0, typename... SliceArgs> struct mask<c,A0,SliceArgs...> { 
     static constexpr bool is_range = std::is_base_of<range, A0>::value; 
     static constexpr ull_t invoke(ull_t m, int P, int c) { return ((is_range ? P+1 : 0) << (4*c)) + mask<SliceArgs...>::invoke(mo,P + is_range, c+1);}
     };
     */

  template <typename... SliceArgs> struct _impl;
  template <>                        struct _impl<> { 
   static constexpr int n_range_ellipsis=0;
   static constexpr ull_t mask(ull_t mo, int P,    int c) { return 0ull;}
   static constexpr ull_t smo (ull_t mo, ull_t ma, int c) { return 0ull;}
  };
  template <typename A0, typename... SliceArgs> struct _impl<A0,SliceArgs...> { 
   static constexpr bool is_range = std::is_base_of<range, A0>::value; // range and ellipsis (derived from range)
   static constexpr int n_range_ellipsis= _impl<SliceArgs...>::n_range_ellipsis + is_range;
   static constexpr ull_t mask(ull_t mo, int P, int c)    { return ((is_range ? P+1ull : 0ull) << (4*(c+1))) + _impl<SliceArgs...>::mask(mo,P + ull_t(is_range), c+1);}
   static constexpr ull_t r2  (ull_t mo, ull_t ma, int c) { return permutations::apply( ma, permutations::apply(permutations::inverse(mo),c)); } 
   static constexpr ull_t aux (ull_t r,  int c)           { return (r==0ull ? ( (r-1ull) << (4*c)) : 0ull);}
   static constexpr ull_t smo (ull_t mo, ull_t ma, int c) { return aux( r2(mo,ma,c), c )  + _impl<SliceArgs...>::smo(mo,ma,c+1); }
  };

  template<typename ... Args> 
   constexpr ull_t sliced_memory_order(ull_t mo) { return _impl<Args...>::smo(mo, _impl<Args...>::mask(mo,0,0),0);}

  template<typename ... Args> 
   constexpr ull_t sliced_memory_order2(ull_t mo) { return 
    (is_c (mo) ? c_order(_impl<Args...>::n_range_ellipsis) : 
     ( is_fortran(mo) ? fortran_order(_impl<Args...>::n_range_ellipsis) : 
       _impl<Args...>::smo(mo, _impl<Args...>::mask(mo,0,0),0) ));}

}}}}//namespace triqs::arrays 
#endif
