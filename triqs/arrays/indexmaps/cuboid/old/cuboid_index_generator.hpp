
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

#ifndef TRIQS_ARRAYS_INDEXMAP_CUBOID_INDEX_GENERATOR_H
#define TRIQS_ARRAYS_INDEXMAP_CUBOID_INDEX_GENERATOR_H
#include "../permutation.hpp"
#include "../../impl/mini_vector.hpp"

namespace triqs { namespace arrays { namespace indexmaps { 

 namespace cuboid_domain_generator_details {  template<typename Perm, int rank, int v> struct index_advance;}
 
 /**
  * Generates the value of the indices of a cuboid_domain.
  * Similar to cuboid_map_iterator, but does not compute a position.
  */
 template <typename CuboidDomain, typename IterationOrder= Permutations::identity<CuboidDomain::rank> >
 class cuboid_index_generator {
   public:
    typedef CuboidDomain domain_type;
    typedef IterationOrder iteration_order;
    typedef typename domain_type::index_value_type indices_type;
    //cuboid_index_generator (): _dom(), atend(true) {}
    cuboid_index_generator (const domain_type & P, bool atEnd=false): _dom(P), atend(atEnd) {}
    cuboid_index_generator & operator++(){ 
     cuboid_domain_generator_details::index_advance<IterationOrder,domain_type::rank,domain_type::rank>::do_it(indices_tuple,_dom,atend);
     return *this;
    }
    
    bool operator==(const cuboid_index_generator & IT2) const { assert((IT2._dom == _dom)); return ((IT2.atend==atend) );}
    //bool operator==(const cuboid_index_generator & IT2) const { return ((IT2._dom == _dom) && (IT2.atend==atend) );}
    bool operator!=(const cuboid_index_generator & IT2) const { return (!operator==(IT2)); }
    indices_type const & operator *() const { return indices_tuple; }
    operator bool () const {return !atend; }
   protected:
    domain_type _dom;
    indices_type indices_tuple; 
    bool atend;
 };

 //  Implementation details 
 namespace cuboid_domain_generator_details {
  template<typename Perm, int rank, int v> struct index_advance {
   static const int p = Permutations::eval<Perm, rank - v >::value;
   template<class I, class D> static void do_it(I & index, D const & _dom, bool & atend) {
#ifdef TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
    if (atend) TRIQS_RUNTIME_ERROR << "Iterator in cuboid can not be pushed after end !";
#endif
    if (index[p] < _dom.lengths()[p]-1) { ++(index[p]); return; }
    index[p] = 0; index_advance<Perm, rank,v-1>::do_it(index, _dom,atend);
   }
  };
  template<typename Perm, int rank> struct index_advance <Perm, rank,0> {
   template<class I, class D> static void do_it(I  & , D const &, bool & atend) { atend = true; }
  };
 }

}}}//namespace triqs::arrays 
#endif
