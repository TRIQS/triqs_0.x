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
#ifndef TRIQS_GF_MATSUBARA_TIME_H
#define TRIQS_GF_MATSUBARA_TIME_H
#include "./tools.hpp"
#include "./gf.hpp"
#include "./local/tail.hpp"
#include "./domains/matsubara.hpp"
#include "./meshes/linear.hpp"

namespace triqs { namespace gf {

 struct imtime {

  /// A tag to recognize the function
  struct tag {};

  /// The domain
  typedef matsubara_domain<false> domain_t;

  /// The Mesh
  typedef linear_mesh<domain_t> mesh_t;

  /// The tail
  typedef local::tail singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

  static std::string h5_name() { return "GfImTime";}
 };

 /// ---------------------------  evaluator ---------------------------------

 template<>
  struct evaluator<imtime> {
   static constexpr int arity = 1;
   //ERROR : give a double and interpolate
   template<typename G>
   arrays::matrix_view<double>  operator() (G const * g,long n)  const {return g->data_view()(n, arrays::range(), arrays::range()); }
   template<typename G>
   local::tail_view operator()(G const * g,freq_infty const &) const {return g->singularity_view();}
  };

 /// ---------------------------  data access  ---------------------------------

 template<> struct data_proxy<imtime> : data_proxy_array<double,3> {};

 // -------------------  ImmutableGfMatsubaraTime identification trait ------------------

 template<typename G> struct ImmutableGfMatsubaraTime : std::is_base_of<typename imtime::tag,G> {};

 // -------------------------------   Factories  --------------------------------------------------

 template<> struct gf_factories< imtime> : imtime { 
  typedef gf<imtime> gf_t;

  static mesh_t make_mesh(double beta, statistic_enum S, size_t n_time_slices, mesh_kind mk) {
   return mesh_t(domain_t(beta,S), 0, beta, n_time_slices, mk);
  }
  static gf_t make_gf(mesh_t && m, tqa::mini_vector<size_t,2> shape, local::tail_view const & t) {
   gf_t::data_non_view_t A(shape.front_append(m.size())); A() =0;
   return gf_t ( m, std::move(A), t, nothing() ) ;
  }
  /*static gf_t make_gf(double beta, statistic_enum S, tqa::mini_vector<size_t,2> shape) {
   return make_gf(make_mesh(beta,S,1025,half_bins), shape, local::tail(shape));
  }
  static gf_t make_gf(double beta, statistic_enum S, tqa::mini_vector<size_t,2> shape, size_t Nmax) {
   return make_gf(make_mesh(beta,S,Nmax,half_bins), shape, local::tail(shape));
  }
  */
  static gf_t make_gf(double beta, statistic_enum S,  tqa::mini_vector<size_t,2> shape, size_t Nmax=1025, mesh_kind mk= half_bins) {
   return make_gf(make_mesh(beta,S,Nmax,mk), shape, local::tail(shape));
  }
  static gf_t make_gf(double beta, statistic_enum S, tqa::mini_vector<size_t,2> shape, size_t Nmax, mesh_kind mk, local::tail_view const & t) {
   return make_gf(make_mesh(beta,S,Nmax,mk), shape, t);
  }
 };
}}
#endif

