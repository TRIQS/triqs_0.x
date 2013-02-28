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
#ifndef TRIQS_GF_ONE_REAL_TIME_H
#define TRIQS_GF_ONE_REAL_TIME_H
#include "./tools.hpp"
#include "./gf.hpp"
#include "./local/tail.hpp"
#include "./meshes/linear.hpp"

namespace triqs { namespace gf { 

 struct retime {

  /// A tag to recognize the function 
  struct tag {};

  /// The domain
  struct domain_t {
   typedef double point_t; 
   bool operator == (domain_t const & D) const { return true; }
   friend void h5_write (h5::group fg, std::string subgroup_name, domain_t const & d) {}
   friend void h5_read  (h5::group fg, std::string subgroup_name, domain_t & d){ }
   friend class boost::serialization::access;
   template<class Archive> void serialize(Archive & ar, const unsigned int version) {}
  };

  /// The Mesh
  typedef linear_mesh<domain_t> mesh_t;
 
  /// The target
  typedef arrays::matrix<std::complex<double> >     target_t;
  //  typedef arrays::matrix<std::complex<double>, arrays::Option::Fortran >     target_t;
  typedef typename target_t::view_type                                       target_view_t;

  /// The storage
  typedef arrays::array<target_t::value_type,3> storage_t;
  typedef typename storage_t::view_type         storage_view_t;

  /// The tail
  typedef local::tail   singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

  /// Indices
  typedef indices_2_t indices_t;

  /// Arity (number of argument in calling the function)
  static const int arity =1;

  struct evaluator { 
   /// All the possible calls of the gf
   template<typename D, typename T>
    target_view_t operator() (mesh_t const & mesh, D const & data, T const & t, double t0)  const {
     size_t index; double w; bool in;
     std::tie(in, index, w) = windowing(mesh,t0);
     if (!in) TRIQS_RUNTIME_ERROR <<" Evaluation out of bounds";
     //return data(arrays::ellipsis(),mesh.index_to_linear(index)); 
     target_t res = w*data(arrays::ellipsis(),mesh.index_to_linear(index)) + (1-w)*data(arrays::ellipsis(),mesh.index_to_linear(index+1));
     return res;
    } 

   template<typename D, typename T>
    local::tail_view operator()(mesh_t const & mesh, D const & data, T const & t, freq_infty const &) const {return t;} 
  };

  struct bracket_evaluator {};

  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS> 
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) { 
    for (size_t u=0; u<mesh.size(); ++u)  { target_view_t( data(tqa::range(),tqa::range(),u)) = rhs(mesh[u]); }
    t = rhs( local::tail::omega(t.shape(),t.size()));
   }

  static std::string h5_name() { return "retime_gf";}

  // -------------------------------   Factories  --------------------------------------------------

  typedef gf<retime> gf_t;

  static gf_t make_gf(double tmin, double tmax, size_t n_time_points, tqa::mini_vector<size_t,2> shape) { 
   retime::mesh_t m(retime::domain_t(), tmin, tmax, n_time_points, mesh_t::full_bins);
   gf_t::data_non_view_t A(shape.append(m.size())); A() =0;
   return gf_t (m, std::move(A), local::tail(shape), nothing(), indices_t(shape) ) ;
  }

 };
 
 // A trait to identify objects that have the concept ImmutableGfOneTime
 template<typename G> struct ImmutableGfOneRealTime : boost::is_base_of<typename retime::tag,G> {};  

}}
#endif

