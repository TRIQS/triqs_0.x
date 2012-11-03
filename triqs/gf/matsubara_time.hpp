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
#include "./gf_proto.hpp"
#include "./domains/matsubara.hpp"
#include "./meshes/linear.hpp"

namespace triqs { namespace gf { 

 struct matsubara_time { 

  /// A tag to recognize the function 
  struct tag {};

  /// The domain
  typedef matsubara_domain<false> domain_t;

  /// The Mesh
  typedef linear_mesh<domain_t> mesh_t;

  /// The target
  typedef arrays::matrix<std::complex<double> >     target_t;
  //  typedef arrays::matrix<std::complex<double>, arrays::Option::Fortran >     target_t;
  typedef typename target_t::view_type            target_view_t;

  /// The tail
  typedef local::tail singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

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
    for (size_t u=0; u<mesh.size(); ++u)  { target_view_t( data(tqa::range(),tqa::range(),u)) = rhs(mesh[u]); }
    t = rhs( local::tail::omega(t.shape(),t.size()));
    // if f is an expression, replace the placeholder with a simple tail. If f is a function callable on freq_infty, 
    // it uses the fact that tail_non_view_t can be casted into freq_infty 
   }

  static std::string h5_name() { return "matsubara_time";}

  // -------------------------------   Factories  --------------------------------------------------

  typedef gf<matsubara_time> gf_t;
  typedef gf_view<matsubara_time> gf_view_t;

  static mesh_t make_mesh (double beta, statistic_enum S, size_t N_time_slices = 1025) {
   return mesh_t( domain_t(beta,S), 0, beta, N_time_slices);
  }

  static gf_t make_gf(mesh_t && m, tqa::mini_vector<size_t,2> shape, local::tail_view const & t) { 
   gf_t::data_non_view_t A(shape.append(m.size())); A() =0;
   return gf_t ( m, std::move(A), t, nothing() ) ;
  }

  static gf_t make_gf(double beta, statistic_enum S, tqa::mini_vector<size_t,2> shape) { 
   return make_gf(make_mesh(beta,S), shape, local::tail(shape));
  }

  static gf_t make_gf(double beta, statistic_enum S, size_t N_time_slices, tqa::mini_vector<size_t,2> shape, local::tail_view const & t) { 
   return make_gf(make_mesh(beta,S,N_time_slices), shape, t);
  }

 };

 // -------------------------------   Expression template --------------------------------------------------

 // A trait to identify objects that have the concept ImmutableGfMatsubaraFreq
 template<typename G> struct ImmutableGfMatsubaraTime : boost::is_base_of<typename matsubara_time::tag,G> {};  

 // This defines the expression template with boost::proto (cf gf_proto.hpp).
 TRIQS_GF_DEFINE_OPERATORS(matsubara_time,local::is_scalar_or_element,ImmutableGfMatsubaraTime);

}}

#endif

