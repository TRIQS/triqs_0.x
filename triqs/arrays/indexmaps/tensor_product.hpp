
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

#ifndef TRIQS_ARRAYS_INDEXMAP_TENSOR_PRODUCT_H
#define TRIQS_ARRAYS_INDEXMAP_TENSOR_PRODUCT_H

namespace triqs { namespace arrays { namespace indexmaps { 

 template <int N, int P> cuboid_domain<N+P> tensor_product( cuboid_domain<N> const & A, cuboid_domain<P> const & B); 

 // the cuboid
 namespace result_of { 
  template<typename T1, typename T2> struct tensor_product; // move it up.
  template <int N, int P> struct tensor_product<cuboid_domain<N>, cuboid_domain<P> > { typedef cuboid_domain<N+P> type;};
 };

 template <int N, int P> 
  cuboid_domain<N+P> tensor_product( cuboid_domain<N> const & A, cuboid_domain<P> const & B) { 
   cuboid_domain<N+P> res; 
   static_assert(0,"not implemented"); // concatenate the tuples.
  }

}}}//namespace triqs::arrays 
#endif

