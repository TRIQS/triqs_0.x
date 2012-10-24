/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_GF_TWO_TIMES_H
#define TRIQS_GF_TWO_TIMES_H

#include <triqs/gf/descriptors/two_times.hpp> 
#include <triqs/gf/gf.hpp>

namespace triqs { namespace gf { 

 /// Make
 gf<two_times> make_gf(two_times, double tmax, double n_time_slices, tqa::mini_vector<size_t,2> shape) { 
  //two_times::mesh_t m(tmax,n_time_slices);
  one_time::mesh_t m1(0, tmax,n_time_slices);
  two_times::mesh_t m(m1,m1);
  gf<two_times>::data_non_view_t A(shape.append(m.size())); A() =0;
  return gf<two_times> (m, std::move(A), nothing(), nothing() ) ;
 }

}}
#endif

