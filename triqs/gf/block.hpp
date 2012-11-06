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
#ifndef TRIQS_GF_BLOCK_H
#define TRIQS_GF_BLOCK_H
#include "./tools.hpp"
#include "./gf.hpp"
#include "./local/tail.hpp"
#include "./gf_proto.hpp"
#include "./meshes/discrete.hpp"

namespace triqs { namespace gf { 

 template<typename Target>
 struct block { 

  /// A tag to recognize the function 
  struct tag {};

  /// The Mesh
  typedef discrete_mesh mesh_t;

  /// The target
  typedef gf<Target> target_t;
  typedef typename target_t::view_type            target_view_t;

  /// 
  typedef nothing singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

  /// Arity (number of argument in calling the function)
  static const int arity =1;

  struct evaluator { 
   template<typename D, typename T>
    target_view_t operator() (mesh_t const & mesh, D const & data, T const & t, long  n)  const {return data[n]; }
  };

  struct bracket_evaluator {
   template<typename D, typename T>
    target_view_t  operator() (mesh_t const & mesh, D & data, T & t, long  n)  const {return data[n]; }
  };

  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS> 
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) { for (auto w: mesh) {data[w.index] = rhs(w); } }

  static std::string h5_name() { return "GFBlock";}

  // -------------------------------   Factories  --------------------------------------------------

  typedef gf<block> gf_t;
  typedef gf_view<block> gf_view_t;

  static gf_t make_gf(std::vector<gf<Target>> const & V)  { return gf_t ( mesh_t(V.size()), V,            nothing(), nothing() ) ; }
  static gf_t make_gf(std::vector<gf<Target>> && V)       { return gf_t ( mesh_t(V.size()), std::move(V), nothing(), nothing() ) ; }

  static gf_view_t make_gf_view(std::vector<gf<Target>> const & V)      { return gf_view_t ( mesh_t(V.size()), V,            nothing(), nothing() ) ; }
  static gf_view_t make_gf_view(std::vector<gf<Target>> && V)           { return gf_view_t ( mesh_t(V.size()), std::move(V), nothing(), nothing() ) ; }
  static gf_view_t make_gf_view(std::vector<gf_view<Target>> const & V) { return gf_view_t ( mesh_t(V.size()), V,            nothing(), nothing() ) ; }
  static gf_view_t make_gf_view(std::vector<gf_view<Target>> && V)      { return gf_view_t ( mesh_t(V.size()), std::move(V), nothing(), nothing() ) ; }

 };


 // -------------------------------   Expression template --------------------------------------------------

 // A trait to identify objects that have the concept ImmutableGfMatsubaraFreq
 //template<typename G> struct ImmutableGfMatsubaraFreq : boost::is_base_of<typename block::tag,G> {};  

 // This defines the expression template with boost::proto (cf gf_proto.hpp).
 //TRIQS_GF_DEFINE_OPERATORS(block,local::is_scalar_or_element,ImmutableGfMatsubaraFreq);

}}

#endif


