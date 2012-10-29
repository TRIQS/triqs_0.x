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
#ifndef TRIQS_GF_DESCRIPTOR_K_MATSUB_H
#define TRIQS_GF_DESCRIPTOR_K_MATSUB_H
#include "../tools.hpp"
#include "../gf.hpp"
#include "../gf_proto.hpp"
#include <array>
#include "./matsubara_freq.hpp"
#include "./BZ.hpp"

namespace triqs { namespace gf { 

 template<typename Omega_t>
 struct k_matsub { 

  /// A tag to recognize the function 
  struct tag {};

  /// The domain
  struct domain_t {
   typedef std::array<double,2> point_t; // a fixed size (2) array of double... name is unfortunate, everything is called array !
   double t_max; 
   domain_t (double t_max_) : t_max(t_max_) {} 
   bool operator == (domain_t const & D) const { return ((std::abs(t_max - D.t_max)<1.e-15) );}
  };

  typedef mesh_product_2<BZ::mesh_t, matsubara_freq::mesh_t> mesh_t;

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

  static typename mesh_t::index_t point_to_index( mesh_t const & mesh, BZ::mesh_t::point_t const & k,  matsubara_freq::mesh_t::point_t const & om) { 
   auto & mesh_k  = mesh.component(zero);
   auto & mesh_om = mesh.component(one);
   return typename mesh_t::index_t( point_to_index(mesh_k,k), point_to_index(mesh_om,om));
  }

static typename mesh_t::index_t point_to_index( mesh_t const & mesh, BZ::mesh_t::point_t const & k,  matsubara_freq::mesh_t::point_t const & om) { 
   return typename mesh_t::index_t( mesh.component(zero).point_to_index(k),  mesh.component(one).point_to_index(om));
  }


  /// All the possible calls of the gf
  template<typename D>
   target_view_t operator() (mesh_t const & mesh, D const & data, nothing , BZ::mesh_t::point_t const & k,  matsubara_freq::mesh_t::point_t const & om)  const {
  //  return data(arrays::range(), arrays::range(),mesh.index_to_linear( point_to_index(mesh,k,om)));
    //return data(arrays::range(), arrays::range(),mesh.index_to_linear( mesh.component(zero).point_to_index(k),  mesh.component(one).point_to_index(om)));
    //
    // renvoie un (i,w) : tuple<int,double>
    int ik; double wk; 
    std::tie(ik,wk) = point_to_index(...);

    return data(arrays::range(), arrays::range(),mesh.index_to_linear( point_to_index(mesh.component(zero), k),  point_to_index(mesh.component(one),om)));
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

static std::string h5_name() { return "k_matsub";}
};

// -------------------------------   Expression template --------------------------------------------------

// A trait to identify objects that have the concept ImmutableGfTwoTimes
template<typename G> struct ImmutableGfK_Matsub : boost::is_base_of<typename k_matsub::tag,G> {};  

// This defines the expression template with boost::proto (cf gf_proto.hpp).
// TRIQS_GF_DEFINE_OPERATORS(two_times,local::is_scalar_or_element,ImmutableGfTwoTimes);

}}
#endif

