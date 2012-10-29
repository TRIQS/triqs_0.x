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
#ifndef TRIQS_GF_DESCRIPTOR_2TIMES_H
#define TRIQS_GF_DESCRIPTOR_2TIMES_H
#include "../tools.hpp"
#include "../gf.hpp"
#include "../gf_proto.hpp"
#include <array>
#include "./one_time.hpp"
#include "../meshes/product.hpp"

namespace triqs { namespace gf { 

 struct two_times { 

  /// A tag to recognize the function 
  struct tag {};

  typedef mesh_product<one_time::mesh_t, one_time::mesh_t> mesh_t;
 
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

  /// Arity (number of argument in calling the function)
  static const int arity =2;

  /// All the possible calls of the gf
  template<typename D, typename T>
   target_view_t operator() (mesh_t const & mesh, D const & data, T const & t, double t0, double t1)  const {
    auto &  m0 = mesh.component(zero); //std::get<0>(mesh.m_tuple);
    //auto &  m0 = std::get<0>(mesh.m_tuple);
    double s= m0.domain().t_max/m0.size();
    return data(arrays::range(), arrays::range(),mesh.index_to_linear( mesh_t::index_t(t0*s, t1*s)));//mesh.index_to_linear(mesh.point_to_index (t1,t2)));
    //return data(arrays::range(), arrays::range(),mesh.index_to_linear( mesh_t::index_t(t0*s, t1*s)));//mesh.index_to_linear(mesh.point_to_index (t1,t2)));
   } 

  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS> 
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) {
    // access to the data . Beware, we view it as a *matrix* NOT an array... (crucial for assignment to scalars !) 
    for (auto & p : mesh) {  
     //std::cout  << " rhs "<< rhs(p.t0,p.t1)<< std::endl;
     //std::cout  << " p "<< p.index[0] << "  "<< p.index[1]<< " linear = "<< p.m->index_to_linear(p.index)<< std::endl;
     target_view_t( data(tqa::range(),tqa::range(),p.m->index_to_linear(p.index))) = rhs(p[zero],p[one]); }
     //target_view_t( data(tqa::range(),tqa::range(),p.m->index_to_linear(p.index))) = rhs(p[ZERO],p[ONE]); }
     //target_view_t( data(tqa::range(),tqa::range(),p.m->index_to_linear(p.index))) = rhs(p[zero_t()],p[one_t()]); }
     //target_view_t( data(tqa::range(),tqa::range(),p.m->index_to_linear(p.index))) = rhs(p._0(),p._1()); }
    //for (size_t u=0; u<mesh.size(); ++u)  { target_view_t( data(tqa::range(),tqa::range(),u)) = rhs(mesh[u]); }
   }

  static std::string h5_name() { return "two_times";}
 };

 // -------------------------------   Expression template --------------------------------------------------

 // A trait to identify objects that have the concept ImmutableGfTwoTimes
 template<typename G> struct ImmutableGfTwoTimes : boost::is_base_of<typename two_times::tag,G> {};  

 // This defines the expression template with boost::proto (cf gf_proto.hpp).
 // TRIQS_GF_DEFINE_OPERATORS(two_times,local::is_scalar_or_element,ImmutableGfTwoTimes);

}}
#endif

