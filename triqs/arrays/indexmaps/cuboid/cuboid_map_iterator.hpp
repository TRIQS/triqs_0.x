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
#ifndef TRIQS_ARRAYS_INDEXMAP_CUBOID_MAP_ITERATOR_H
#define TRIQS_ARRAYS_INDEXMAP_CUBOID_MAP_ITERATOR_H
#include <boost/iterator/iterator_facade.hpp>
#include "../../impl/mini_vector.hpp"

namespace triqs { namespace arrays { namespace indexmaps { 

 /**
  * Iterator on a cuboid_map, modeling the IndexMapIterator concept.
  * Iteration order is the order in which to iterate on the indices.
  * It is given by a permutation, with the same convention as IndexOrder.
  */
 template <typename CuboidMap, typename IterationOrder= Permutations::identity<CuboidMap::rank> >
  class cuboid_map_iterator :  
   public boost::iterator_facade< cuboid_map_iterator<CuboidMap,IterationOrder> , const std::ptrdiff_t , boost::forward_traversal_tag > {
    public:
     typedef CuboidMap indexmap_type;
     typedef typename indexmap_type::domain_type domain_type;
     typedef IterationOrder iteration_order;
     typedef typename domain_type::index_value_type indices_type;
     typedef const std::ptrdiff_t return_type;

     //cuboid_map_iterator (): Parent(NULL), pos(0),atend(true) {}
     cuboid_map_iterator (const indexmap_type & P, bool atEnd=false): im(P), pos(im.start_shift()),atend(atEnd) {}

     indices_type const & indices() const { return indices_tuple; }
     operator bool() const { return !atend;}

    private:
     friend class boost::iterator_core_access;
     void increment(); // below
     // I assume that the two iterator are on the same cuboid_map
     bool equal(cuboid_map_iterator const & other) const {return ((other.atend==atend)&&(other.pos==pos));}
     //bool equal(cuboid_map_iterator const & other) const {return ((other.im==im)&&(other.atend==atend)&&(other.pos==pos));}
     return_type & dereference() const { assert (!atend); return pos; }

     indexmap_type im;
     indices_type indices_tuple; 
     std::ptrdiff_t pos;
     bool atend;
   };

 /*-----------------------------------------
  *       Implementation details 
  *----------------------------------------*/

 namespace cuboid_map_iterator_details {
  template<typename CM, typename IO, int rank, int v> struct inc_ind_impl {
   static const int p = Permutations::eval<IO, rank - v >::value;
   typedef mini_vector<size_t,CM::rank> indices_type;
   static void do_it(indices_type & index, CM const & im, std::ptrdiff_t & pos, bool & atend) {
#ifdef TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
    if (atend) TRIQS_RUNTIME_ERROR << "Iterator in cuboid can not be pushed after end !";
#endif
    if (index[p] < im.lengths()[p]-1) { ++(index[p]); pos += im.strides()[p]; return; }
    index[p] = 0; 
    pos -= (im.lengths()[p]-1) * im.strides()[p];
    inc_ind_impl<CM,IO, rank,v-1>::do_it(index, im,pos,atend);
   }
  };

  template<typename CM, typename IO, int rank> struct inc_ind_impl <CM,IO, rank,0> {
   typedef mini_vector<size_t,CM::rank> indices_type;
   static void do_it(indices_type & , CM const & ,std::ptrdiff_t & , bool & atend) { atend = true; }
  };
 }

 template <typename CM, typename IO>
  void cuboid_map_iterator<CM,IO>::increment() { 
   cuboid_map_iterator_details::inc_ind_impl<CM,IO,CM::rank,CM::rank>::do_it(indices_tuple,im,pos,atend);
  }

}}}//namespace triqs::arrays 
#endif
