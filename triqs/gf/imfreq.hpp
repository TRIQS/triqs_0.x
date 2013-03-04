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
#ifndef TRIQS_GF_MATSUBARA_FREQ_H
#define TRIQS_GF_MATSUBARA_FREQ_H
#include "./tools.hpp"
#include "./gf.hpp"
#include "./local/tail.hpp"
#include "./domains/matsubara.hpp"
#include "./meshes/linear.hpp"

namespace triqs { namespace gf {

 struct imfreq {

  /// A tag to recognize the function
  struct tag {};

  /// The domain
  typedef matsubara_domain<true> domain_t;

  /// The Mesh
  typedef linear_mesh<domain_t> mesh_t;

  /// The tail
  typedef local::tail singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

  /// Indices
  typedef indices_2_t indices_t;

  static std::string h5_name() { return "imfreq_gf";}

 };

 /// ---------------------------  evaluator ---------------------------------

 template<typename G>
  struct evaluator<imfreq,G> {
   static const int arity =1;/// Arity (number of argument in calling the function)
   G const * g; evaluator(G const & g_): g(&g_){}
   arrays::matrix_view<std::complex<double> >  operator() (long n)  const {return g->data_view()(arrays::range(), arrays::range(),n); }
   local::tail_view operator()(freq_infty const &) const {return g->singularity_view();}
  };

 /// ---------------------------  data access  ---------------------------------

 template<> struct data_proxy<imfreq> : data_proxy_array<std::complex<double>,3> {};

 // -------------------  ImmutableGfMatsubaraFreq identification trait ------------------
 
 template<typename G> struct ImmutableGfMatsubaraFreq : boost::is_base_of<typename imfreq::tag,G> {};

 // -------------------------------   Factories  --------------------------------------------------

 template<> struct gf_factories<imfreq> : imfreq { 
  typedef gf<imfreq> gf_t;

  static mesh_t make_mesh (double beta, statistic_enum S, size_t Nmax = 1025) {
   double m1 = std::acos(-1)/beta;
   return mesh_t( domain_t(beta,S), m1, (2*Nmax+1)*m1, Nmax, without_last);
  }

  static gf_t make_gf(mesh_t && m, tqa::mini_vector<size_t,2> shape, local::tail_view const & t) {
   gf_t::data_non_view_t A(shape.append(m.size())); A() =0;
   //gf_t::data_non_view_t A(shape.append(m.size()), FORTRAN_LAYOUT); A() =0;
   return gf_t ( m, std::move(A), t, nothing(), indices_t(shape) ) ;
  }

  static gf_t make_gf(double beta, statistic_enum S, tqa::mini_vector<size_t,2> shape) {
   return make_gf(make_mesh(beta,S), shape, local::tail(shape));
  }

  static gf_t make_gf(double beta, statistic_enum S,  tqa::mini_vector<size_t,2> shape, size_t Nmax) {
   return make_gf(make_mesh(beta,S,Nmax), shape, local::tail(shape));
  }

  static gf_t make_gf(double beta, statistic_enum S, tqa::mini_vector<size_t,2> shape, size_t Nmax, local::tail_view const & t) {
   return make_gf(make_mesh(beta,S,Nmax), shape, t);
  }

 };
}}

#endif

