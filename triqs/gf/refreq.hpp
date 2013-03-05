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
#ifndef TRIQS_GF_FREQ_H
#define TRIQS_GF_FREQ_H
#include "./tools.hpp"
#include "./gf.hpp"
#include "./local/tail.hpp"
#include "./domains/R.hpp"
#include "./meshes/linear.hpp"

// Shall we use the same type as for retime : same code, almost ??

namespace triqs { namespace gf {

 struct refreq {

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

  static std::string h5_name() { return "refreq_gf";}
 };

 /// ---------------------------  evaluator ---------------------------------

 template<>
  struct evaluator<refreq> {
   template<typename G>
    arrays::matrix_view<std::complex<double> >  operator() (G const * g,double w0)  const {
     auto & data = g->data_view();
     auto & mesh = g->mesh();
     size_t index; double w; bool in;
     std::tie(in, index, w) = windowing(mesh,w0);
     if (!in) TRIQS_RUNTIME_ERROR <<" Evaluation out of bounds";
     arrays::matrix<std::complex<double> > res = w*data(arrays::ellipsis(),mesh.index_to_linear(index)) + (1-w)*data(arrays::ellipsis(),mesh.index_to_linear(index+1));
     return res;
    }
   template<typename G>
    local::tail_view operator()(G const * g,freq_infty const &) const {return g->singularity_view();}
  };

 /// ---------------------------  data access  ---------------------------------

 template<> struct data_proxy<refreq> : data_proxy_array<std::complex<double>,3> {};

 // -------------------  ImmutableGfFreq identification trait ------------------

 template<typename G> struct ImmutableGfFreq : boost::is_base_of<typename refreq::tag,G> {};

 // -------------------------------   Factories  --------------------------------------------------

 template<> struct gf_factories<refreq> : refreq { 
  typedef gf<refreq> gf_t;

  static mesh_t make_mesh(double wmin, double wmax, size_t n_freq, mesh_kind mk) {
   return mesh_t(domain_t(), wmin, wmax, n_freq, mk);
  }

  static gf_t make_gf(double wmin, double wmax, size_t n_freq, tqa::mini_vector<size_t,2> shape) {
   gf_t::data_non_view_t A(shape.append(n_freq)); A() =0;
   return gf_t(make_mesh(wmin, wmax, n_freq, full_bins), std::move(A), local::tail(shape), nothing(), indices_t(shape));
  }

  static gf_t make_gf(double wmin, double wmax, size_t n_freq, tqa::mini_vector<size_t,2> shape, mesh_kind mk) {
   gf_t::data_non_view_t A(shape.append(n_freq)); A() =0;
   return gf_t(make_mesh(wmin, wmax, n_freq, mk), std::move(A), local::tail(shape), nothing(), indices_t(shape));
  }

 };

}}
#endif

