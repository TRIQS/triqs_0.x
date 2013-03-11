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
#include "./domains/R.hpp"
#include "./meshes/linear.hpp"

namespace triqs { namespace gf {

 struct retime {

  /// A tag to recognize the function
  struct tag {};

  /// The domain
  typedef R_domain domain_t;

  /// The Mesh
  typedef linear_mesh<domain_t> mesh_t;

  /// The tail
  typedef local::tail   singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

  /// Indices
  typedef indices_2_t indices_t;

  static std::string h5_name() { return "retime_gf";}

 };

 /// ---------------------------  evaluator ---------------------------------

 template<>
  struct evaluator<retime> {
   template<typename G>
    arrays::matrix_view<std::complex<double> >  operator() (G const * g,double t0)  const {
     auto & data = g->data_view();
     auto & mesh = g->mesh();
     size_t index; double w; bool in;
     std::tie(in, index, w) = windowing(mesh,t0);
     if (!in) TRIQS_RUNTIME_ERROR <<" Evaluation out of bounds";
     arrays::matrix<std::complex<double> > res = w*data(arrays::ellipsis(),mesh.index_to_linear(index)) + (1-w)*data(arrays::ellipsis(),mesh.index_to_linear(index+1));
     return res;
    }
   template<typename G>
    local::tail_view operator()(G const * g,freq_infty const &) const {return g->singularity_view();}
  };

 /// ---------------------------  data access  ---------------------------------

 template<> struct data_proxy<retime> : data_proxy_array<std::complex<double>,3> {};

 // -------------------  ImmutableGfOneRealTime identification trait ------------------

 template<typename G> struct ImmutableGfOneRealTime : boost::is_base_of<typename retime::tag,G> {};

 // -------------------------------   Factories  --------------------------------------------------

 template<> struct gf_factories<retime> : retime { 
  typedef gf<retime> gf_t;

  static mesh_t make_mesh(double tmin, double tmax, size_t n_time_points, mesh_kind mk) {
   return mesh_t(domain_t(), tmin, tmax, n_time_points, mk);
  }

  static gf_t make_gf(double tmin, double tmax, size_t n_time_points, tqa::mini_vector<size_t,2> shape) {
   gf_t::data_non_view_t A(shape.append(n_time_points)); A() =0;
   return gf_t(make_mesh(tmin, tmax, n_time_points, full_bins), std::move(A), local::tail(shape), nothing(), indices_t(shape));
  }

  static gf_t make_gf(double tmin, double tmax, size_t n_time_points, tqa::mini_vector<size_t,2> shape, mesh_kind mk) {
   gf_t::data_non_view_t A(shape.append(n_time_points)); A() =0;
   return gf_t(make_mesh(tmin, tmax, n_time_points, mk), std::move(A), local::tail(shape), nothing(), indices_t(shape));
  }

 };

}}
#endif

