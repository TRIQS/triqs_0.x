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
#ifndef TRIQS_ARRAYS_INDEXMAP_MEMORY_LAYOUT_H
#define TRIQS_ARRAYS_INDEXMAP_MEMORY_LAYOUT_H
#include "../permutation.hpp"
namespace triqs { namespace arrays {

 namespace indexmaps { namespace mem_layout {
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

 }}
 
 struct memory_layout_fortran {};
 struct memory_layout_c {};
 
#define FORTRAN_LAYOUT (triqs::arrays::memory_layout_fortran())
 //struct custom_memory_layout {};

 // stores the layout == order of the indices in memory
 // wrapped into a little type to make constructor unambigous.
 template<int Rank>
  struct memory_layout { 
   ull_t value; 
   memory_layout() { value =indexmaps::mem_layout::c_order(Rank);} 
   explicit memory_layout (ull_t v) : value(v) {assert((permutations::size(mo)==Rank));} 
   memory_layout (const char * ml) {
    assert( (ml[0]=='C') || (ml[0] == 'F'));
    value = (ml[0]=='F' ? indexmaps::mem_layout::fortran_order(Rank) : indexmaps::mem_layout::c_order(Rank));
    //assert( (ml=='C') || (ml == 'F'));
    //value = (ml=='F' ? indexmaps::mem_layout::fortran_order(Rank) : indexmaps::mem_layout::c_order(Rank));
   }
   memory_layout (memory_layout_fortran) { value = indexmaps::mem_layout::fortran_order(Rank); }
   memory_layout (memory_layout_c) { value = indexmaps::mem_layout::c_order(Rank); }
   template<typename ... INT>
    explicit memory_layout(int i0, int i1, INT ... in) : value (permutations::permutation(i0,i1,in...)){
    static_assert( sizeof...(in)==Rank-2, "Error");
    }
   bool operator ==( memory_layout const & ml) const { return value == ml.value;}
   bool operator !=( memory_layout const & ml) const { return value != ml.value;}
  };

}}//namespace triqs::arrays 
#endif
