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
#include "./gf_proto.hpp"
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

  /// The tail
  typedef nothing singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;
 
  /// Indices
  typedef indices_2_t indices_t;

  /// Arity (number of argument in calling the function)
  static const int arity =1;

  /// All the possible calls of the gf
  //ERROR : give a double and interpolate
  struct evaluator { 
  template<typename D, typename T>
   target_view_t operator() (mesh_t const & mesh, D const & data, T const & t, long  n)  const {return data(arrays::range(), arrays::range(),n); }

  template<typename D, typename T>
   local::tail_view operator()(mesh_t const & mesh, D const & data, T const & t, freq_infty const &) const {return t;} 
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

 // -------------------------------   Expression template --------------------------------------------------

 // A trait to identify objects that have the concept ImmutableGfMatsubaraFreq
 template<typename G> struct ImmutableGfLegendre : boost::is_base_of<typename legendre::tag,G> {};  

 // This defines the expression template with boost::proto (cf gf_proto.hpp).
 TRIQS_GF_DEFINE_OPERATORS(legendre,legendre::tag, 1, local::is_scalar_or_element, ImmutableGfLegendre);

}}

#endif

