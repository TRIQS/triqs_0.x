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
#ifndef TRIQS_GF_MATSUBARA_H
#define TRIQS_GF_MATSUBARA_H

#include <triqs/gf/descriptors/matsubara_freq.hpp> 
#include <triqs/gf/descriptors/matsubara_time.hpp> 

#include <triqs/gf/gf.hpp>

namespace triqs { namespace gf { 

 /// Make
 gf<matsubara_freq> make_gf(matsubara_freq, matsubara_freq::mesh_t const & m, tqa::mini_vector<size_t,2> shape, local::tail_view const & t) { 
  gf<matsubara_freq>::data_non_view_t A(shape.append(m.size())); A() =0;
  return gf<matsubara_freq> ( m, std::move(A), t, nothing() ) ;
 }

 /// Make
 gf<matsubara_freq> make_gf(matsubara_freq, double beta, statistic_enum S, tqa::mini_vector<size_t,2> shape) { 
  return make_gf(matsubara_freq(), matsubara_freq::mesh_t(beta,S), shape,local::tail(shape));
 }
 
 /// Make
 gf<matsubara_time> make_gf(matsubara_time, matsubara_time::mesh_t const & m, tqa::mini_vector<size_t,2> shape, local::tail_view const & t) { 
  gf<matsubara_time>::data_non_view_t A(shape.append(m.size())); A() =0;
  return gf<matsubara_time> ( m, std::move(A), t, nothing() ) ;
 }

 /// Make
 gf<matsubara_time> make_gf(matsubara_time, double beta, statistic_enum S, tqa::mini_vector<size_t,2> shape) { 
  return make_gf(matsubara_time(), matsubara_time::mesh_t(beta,S), shape,local::tail(shape));
 }
 
}}
#endif

