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
#include "./meshes/linear.hpp"

// Shall we use the same type as for retime : same code, almost ??

namespace triqs { namespace gf {

 struct refreq {

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

  /// The storage
  typedef arrays::array<std::complex<double>,3> storage_t;
  typedef typename storage_t::view_type         storage_view_t;

  /// The tail
  typedef local::tail   singularity_t;

  /// Symmetry
  typedef nothing symmetry_t;

  /// Indices
  typedef indices_2_t indices_t;

  /// ---------------------------  evaluator ---------------------------------

  template<typename G>
   struct evaluator {
    static const int arity =1;/// Arity (number of argument in calling the function)
    G const * g; evaluator(G const & g_): g(&g_){}
    arrays::matrix_view<std::complex<double> >  operator() (double w0)  const {
     auto & data = g->data_view();
     auto & mesh = g->mesh();
     size_t index; double w; bool in;
     std::tie(in, index, w) = windowing(mesh,w0);
     if (!in) TRIQS_RUNTIME_ERROR <<" Evaluation out of bounds";
      arrays::matrix<std::complex<double> > res = w*data(arrays::ellipsis(),mesh.index_to_linear(index)) + (1-w)*data(arrays::ellipsis(),mesh.index_to_linear(index+1));
     return res;
    }
    local::tail_view operator()(freq_infty const &) const {return g->singularity_view();}
   };

  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS>
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) {
    for (size_t u=0; u<mesh.size(); ++u)  { target_view_t( data(tqa::range(),tqa::range(),u)) = rhs(mesh[u]); }
    t = rhs( local::tail::omega(t.shape(),t.size()));
   }

  static std::string h5_name() { return "refreq_gf";}

  // -------------------------------   Factories  --------------------------------------------------

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

 // A trait to identify objects that have the concept ImmutableGfFreq
 template<typename G> struct ImmutableGfFreq : boost::is_base_of<typename refreq::tag,G> {};
}}
#endif

