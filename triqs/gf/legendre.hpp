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
#ifndef TRIQS_GF_LEGENDRE_TIME_H
#define TRIQS_GF_LEGENDRE_TIME_H
#include "./tools.hpp"
#include "./gf.hpp"
#include "./local/tail.hpp"
#include "./domains/legendre.hpp"
#include "./meshes/discrete.hpp"

namespace triqs { namespace gf {

 struct legendre {

  /// A tag to recognize the function
  struct tag {};

  /// The domain
  typedef legendre_domain domain_t;

  /// The Mesh
  typedef discrete_mesh<domain_t> mesh_t;

  /// The tail
  typedef nothing singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

  /// Arity (number of argument in calling the function)
  static const int arity =1;

  static std::string h5_name() { return "legendre_gf";}

 };
 /// ---------------------------  evaluator ---------------------------------

 template<>
  struct evaluator<legendre> {
   static constexpr int arity = 1;
   //ERROR : give a double and interpolate
   template<typename G>
   arrays::matrix_view<double >  operator() (G const * g,long n)  const {return g->data_view()(n, arrays::range(), arrays::range()); }
   template<typename G>
   local::tail_view operator()(G const * g,freq_infty const &) const {return g->singularity_view();}
  };

 /// ---------------------------  data access  ---------------------------------

 template<> struct data_proxy<legendre> : data_proxy_array<double,3> {};

 // -------------------  ImmutableGfLegendre identification trait ------------------

 template<typename G> struct ImmutableGfLegendre : boost::is_base_of<typename legendre::tag,G> {};

 // -------------------------------   Factories  --------------------------------------------------

 template<> struct gf_factories< legendre> : legendre { 
  typedef gf<legendre> gf_l;

  static mesh_t make_mesh(double beta, statistic_enum S, size_t n_leg) {
   return mesh_t(domain_t(beta,S,n_leg));
  }

  static gf_l make_gf(double beta, statistic_enum S, size_t n_leg, tqa::mini_vector<size_t,2> shape) {
   gf_l::data_non_view_t A(shape.front_append(n_leg)); A() = 0;
   return gf_l(make_mesh(beta, S, n_leg), std::move(A), nothing(), nothing());
  }

 };

}}
#endif

