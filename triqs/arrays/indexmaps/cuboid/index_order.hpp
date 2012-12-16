/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by O. Parcollet
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
  /* The storage order is given by a permutation P stored in a ull_t (unsigned long long) as in permutations::..
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
   */
//#define TRIQS_WORKAROUND_INTEL_COMPILER_BUGS2
#ifndef TRIQS_WORKAROUND_INTEL_COMPILER_BUGS2

  template <typename... SliceArgs> struct _impl{ 
   static constexpr int n_range_ellipsis=0;
   static constexpr ull_t mask(ull_t mo, int P, int c) { return 0;}
   static constexpr ull_t smo (ull_t mo, ull_t ma, int P, int c) { return 0;}
  };
  template <typename A0, typename... SliceArgs> struct _impl< A0,SliceArgs...> { 
   typedef _impl<SliceArgs...> next_t;
   static constexpr bool is_range = std::is_base_of<range, A0>::value; // range and ellipsis (derived from range)
   static constexpr int n_range_ellipsis= next_t::n_range_ellipsis + is_range;
   static constexpr ull_t mask(ull_t mo, int P=0, int c=0)            { return ((is_range ? P+1 : 0) << (4*c)) + next_t::mask(mo,P + is_range,c+1);}
   static constexpr ull_t r2  (ull_t mo, ull_t ma, int c)             { return (ma >> (4*(permutations::apply(permutations::inverse(mo),c)))) & 0xF; } 
   static constexpr ull_t aux (ull_t mo,ull_t ma,ull_t r,int P,int c) { return (r!=0? ( (r-1) << (4*P)) : 0) + next_t::smo(mo,ma,P+(r!=0), c+1); }
   static constexpr ull_t smo (ull_t mo, ull_t ma, int P=0, int c=0)  { return aux(mo,ma,r2(mo,ma,c),P,c); }
  };

  template<typename ... Args> 
   constexpr ull_t sliced_memory_order1(ull_t mo) { return _impl<Args...>::n_range_ellipsis + 0x10*_impl<Args...>::smo(mo, _impl<Args...>::mask(mo));}

  template<typename ... Args> 
   constexpr ull_t sliced_memory_order2(ull_t mo) { return (is_c (mo) ? c_order(_impl<Args...>::n_range_ellipsis) : 
     ( is_fortran(mo) ? fortran_order(_impl<Args...>::n_range_ellipsis) : sliced_memory_order1<Args...>(mo) ));
   }

  template<ull_t mo, typename ... Args> struct sliced_memory_order { static constexpr ull_t value = sliced_memory_order1<Args...>(mo); };

#else
  template <int P, int c, typename... SliceArgs> struct mask { 
   static constexpr ull_t value = 0;
   static constexpr int n_range_ellipsis= 0;
  };
  template <int P, int c, typename A0, typename... SliceArgs> struct mask<P,c,A0,SliceArgs...> { 
   static constexpr bool is_range = std::is_base_of<range, A0>::value; 
   typedef mask<P + is_range, c+1,SliceArgs...> next_t;
   static constexpr int n_range_ellipsis= next_t::n_range_ellipsis + is_range;
   static constexpr ull_t value = ((is_range ? P+1 : 0) << (4*c)) + next_t::value;
  };

  template <ull_t mo, ull_t ma, int P, int c, typename ... A> struct _impl { static constexpr ull_t value = 0; };
  template <ull_t mo, ull_t ma, int P, int c, typename A0, typename... SliceArgs> struct _impl<mo,ma,P,c, A0,SliceArgs...> { 
   static constexpr bool is_range = std::is_base_of<range, A0>::value; 
   struct r2 { // Intel .... 
    static constexpr ull_t pn = permutations::apply(permutations::inverse(mo),c);
    static constexpr ull_t value  = (ma >> (4*(pn))) & 0xFull;
   };
   static const bool r2_nozero = (r2::value!=0);
   static constexpr ull_t value =  (r2::value!=0  ? ( (r2::value-1) << (4*P) ) : 0)+ _impl<mo,ma,P+r2_nozero,c+1,SliceArgs...>::value;
  };

  template<ull_t mo, typename ... Args> struct sliced_memory_order { 
   typedef mask<0,0,Args...> mask_t;
   static constexpr ull_t value = mask_t::n_range_ellipsis + 0x10* _impl<mo, mask_t::value,0,0,Args...>::value ;
  };
#endif
}}}}//namespace triqs::arrays 
#endif
