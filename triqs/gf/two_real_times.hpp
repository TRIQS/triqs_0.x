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
#include "./tools.hpp"
#include "./gf.hpp"
#include "./gf_proto.hpp"
#include "./retime.hpp"
#include "./meshes/product.hpp"

namespace triqs { namespace gf { 

 struct two_times { 

  /// A tag to recognize the function 
  struct tag {};

  typedef mesh_product<retime::mesh_t, retime::mesh_t> mesh_t;

  // suppress from the concept : can always be deduced ? 
  typedef typename mesh_t::domain_t domain_t;

  /// The target
  typedef arrays::matrix<std::complex<double> >     target_t;
  //  typedef arrays::matrix<std::complex<double>, arrays::Option::Fortran >     target_t;
  typedef typename target_t::view_type                                       target_view_t;

  /// 
  typedef nothing singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

  /// Indices
  typedef nothing indices_t;

  /// Arity (number of argument in calling the function)
  static const int arity =2;

  /// All the possible calls of the gf
  struct evaluator { 
  template<typename D, typename T>
   target_view_t operator() (mesh_t const & mesh, D const & data, T const & t, double t0, double t1)  const {
    auto & m0 = mesh.component(zero); //std::get<0>(mesh.m_tuple);
    //auto &  m0 = std::get<0>(mesh.m_tuple);
    double s= m0.x_max()/m0.size();
    return data(arrays::range(), arrays::range(),mesh.index_to_linear( mesh_t::index_t(t0*s, t1*s)));//mesh.index_to_linear(mesh.point_to_index (t1,t2)));
   } 
  };

  struct bracket_evaluator {};
  
  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS> 
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) {
    // access to the data . Beware, we view it as a *matrix* NOT an array... (crucial for assignment to scalars !) 
    for (auto & p : mesh) {  
     target_view_t( data(tqa::range(),tqa::range(),p.m->index_to_linear(p.index))) = rhs(p[zero],p[one]); 
    }
   }

  static std::string h5_name() { return "two_times";}

  // -------------------------------   Factories  --------------------------------------------------

  typedef gf<two_times> gf_t;

  static gf_t make_gf(double tmax, double n_time_slices, tqa::mini_vector<size_t,2> shape) { 
   retime::mesh_t m1(retime::domain_t(),0, tmax,n_time_slices);
   mesh_t m(m1,m1);
   gf_t::data_non_view_t A(shape.append(m.size())); A() =0;
   return gf_t (m, std::move(A), nothing(), nothing() ) ;
  }

 };

 // -------------------------------   Expression template --------------------------------------------------

 // A trait to identify objects that have the concept ImmutableGfTwoTimes
 template<typename G> struct ImmutableGfTwoTimes : boost::is_base_of<typename two_times::tag,G> {};  

 // This defines the expression template with boost::proto (cf gf_proto.hpp).
 // TRIQS_GF_DEFINE_OPERATORS(two_times,2,local::is_scalar_or_element,ImmutableGfTwoTimes);

}}
#endif

