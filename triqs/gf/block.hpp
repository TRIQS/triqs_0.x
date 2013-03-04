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
#include "./meshes/discrete.hpp"

namespace triqs { namespace gf {

 struct block_tag {};

 template<typename Target>
  struct block {

   /// A tag to recognize the function
   struct tag : block_tag {};

   /// The Mesh
   typedef discrete_mesh<discrete_domain> mesh_t;

   /// The storage
   typedef std::vector<gf<Target>> storage_t;
   typedef std::vector<gf_view<Target>> storage_view_t;

   ///
   typedef nothing singularity_t;

   /// Symmetry
   typedef nothing symmetry_t;

   /// Indices
   typedef nothing indices_t;

 
  /// How to fill a gf from an expression (RHS)
  template<typename D, typename T, typename RHS>
   static void assign_from_expression (mesh_t const & mesh, D & data, T & t, RHS rhs) { for (auto w: mesh) {data[w.index] = rhs(w); } }

  static std::string h5_name() { return "block_gf";}

  typedef void has_special_h5_read_write_tag;

  static void h5_data_write(h5::group g, std::string const & s, mesh_t const & mesh, storage_t const & data) {
   auto gr =  g.create_group(s);
   for (size_t i =0; i<mesh.size(); ++i) h5_write(gr,mesh.domain().names()[i],data[i]);
  }

  static void h5_data_read(h5::group g, std::string const & s, mesh_t const & mesh, storage_t & data) {
   auto gr =  g.create_group(s);
   for (size_t i =0; i<mesh.size(); ++i) h5_write(gr,mesh.domain().names()[i],data[i]);
  }
  };
  /// ---------------------------  evaluator ---------------------------------

  template<typename Target, typename G>  
   struct evaluator<block<Target>,G> { 
    static const int arity =1;/// Arity (number of argument in calling the function)
    G const * g; evaluator(G const & g_): g(&g_){}
    //gf_view<Target> operator() (long  n)  const {return g->data_view()[n]; }
   };


  // -------------------------------   Factories  --------------------------------------------------

 template<typename Target>
  struct gf_factories<block<Target>> : block<Target> { 
   typedef block<Target> B;
   typedef typename B::mesh_t mesh_t;
   typedef gf<block<Target>> gf_t;
   typedef gf_view<block<Target>> gf_view_t;

   static gf_t make_gf(std::vector<gf<Target>> const & V)  { return gf_t ( mesh_t(V.size()), V,            nothing(), nothing() ) ; }
   static gf_t make_gf(std::vector<gf<Target>> && V)       { return gf_t ( mesh_t(V.size()), std::move(V), nothing(), nothing() ) ; }

   static gf_t make_gf(std::vector<std::string> const & block_names, std::vector<gf<Target>> const & V) {
    return gf_t(mesh_t(block_names), V, nothing(), nothing() );
   }
   static gf_t make_gf(std::vector<std::string> const & block_names, std::vector<gf<Target>> && V) {
    return gf_t(mesh_t(block_names), std::move(V), nothing(), nothing() );
   }

   /* template<typename... Args>
      static gf_t make_gf(size_t N, Args&& ...args)  {
      std::vector<gf<Target>> V; V.reserve(N);
      for (size_t i=0; i<N; ++i) V.push_back( Target::make_gf (std::forward<Args>(args...)));
      return make_gf(V);
      }
      */
   static gf_t make_gf(int N, gf<Target> const & g)  {
    std::vector<gf<Target>> V; V.reserve(N);
    for (size_t i=0; i<N; ++i)  V.push_back(g);
    return make_gf(V);
   }

   static gf_t make_gf(std::vector<std::string> const & block_names, gf<Target> const & g)  {
    std::vector<gf<Target>> V; V.reserve(block_names.size());
    for (size_t i=0; i<block_names.size(); ++i)  V.push_back(g);
    return make_gf(block_names,V);
   }


   /*  template<typename... Args>
       static gf_t make_gf(std::vector<std::string> const & block_names, Args&& ...args)  {
       std::vector<gf<Target>> V; V.reserve(block_names.size());
       for (size_t i=0; i<block_names.size(); ++i)  V.push_back( Target::make_gf (std::forward<Args>(args...)));
       return make_gf(block_names,V);
       }
       */

   template<typename GF>
    static gf_view_t make_gf_view(std::vector<GF> const & V) { return gf_view_t ( mesh_t(V.size()), V,            nothing(), nothing() ) ; }
   template<typename GF>
    static gf_view_t make_gf_view(std::vector<GF> && V)      { return gf_view_t ( mesh_t(V.size()), std::move(V), nothing(), nothing() ) ; }

  };


  // A trait to identify objects that have the concept ImmutableGfMatsubaraFreq
  template<typename G> struct ImmutableBlockGf : boost::is_base_of<block_tag,G> {};

}}

#endif


