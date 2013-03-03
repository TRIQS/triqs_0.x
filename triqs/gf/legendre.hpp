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

  /// The target
  typedef arrays::matrix<double>         target_t;
  typedef typename target_t::view_type   target_view_t;

  /// The storage
  typedef arrays::array<target_t::value_type,3> storage_t;
  typedef typename storage_t::view_type         storage_view_t;

  /// The tail
  typedef nothing singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

  /// Indices
  typedef indices_2_t indices_t;

  /// Arity (number of argument in calling the function)
  static const int arity =1;

 /// ---------------------------  evaluator ---------------------------------

  template<typename G>
   struct evaluator {
    static const int arity =1;/// Arity (number of argument in calling the function)
    //ERROR : give a double and interpolate
    G const * g; evaluator(G const & g_): g(&g_){}
    arrays::matrix_view<double >  operator() (long n)  const {return g->data_view()(arrays::range(), arrays::range(),n); }
    local::tail_view operator()(freq_infty const &) const {return g->singularity_view();}
   };

  struct bracket_evaluator {};

  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS>
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) {
    // access to the data . Beware, we view it as a *matrix* NOT an array... (crucial for assignment to scalars !)
    for (size_t u=0; u<mesh.size(); ++u)  { target_view_t(data(tqa::range(),tqa::range(),u)) = rhs(mesh[u]); }
   }

  static std::string h5_name() { return "legendre_gf";}

  // -------------------------------   Factories  --------------------------------------------------

  typedef gf<legendre> gf_l;

  static mesh_t make_mesh(double beta, statistic_enum S, size_t n_leg) {
   return mesh_t(domain_t(beta,S,n_leg));
  }

  static gf_l make_gf(double beta, statistic_enum S, size_t n_leg, tqa::mini_vector<size_t,2> shape) {
   gf_l::data_non_view_t A(shape.append(n_leg)); A() = 0;
   return gf_l(make_mesh(beta, S, n_leg), std::move(A), nothing(), nothing(), indices_t(shape));
  }

 };

 // A trait to identify objects that have the concept ImmutableGfMatsubaraFreq
 template<typename G> struct ImmutableGfLegendre : boost::is_base_of<typename legendre::tag,G> {};
}}

#endif

