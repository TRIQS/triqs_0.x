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
#ifndef TRIQS_GF_MESH_LINEAR_H
#define TRIQS_GF_MESH_LINEAR_H
#include "../tools.hpp"
namespace triqs { namespace gf { 

 template<typename Domain>
  struct linear_mesh {

   typedef Domain domain_t;
   typedef size_t index_t; 
   typedef typename domain_t::point_t  domain_pt_t;

   mesh_t (domain_t && dom, domain_pt_t x_min, domain_pt_t x_max, size_t ns) : _dom(dom),L(ns), xmin(x_min), xmax(x_max), step( (x_max - x_min) /(ns-1)){}
   mesh_t () : _dom(),L(0){}

   domain_t const & domain() const { return _dom;}
   size_t size() const {return L;}

   /// Conversions point <-> index <-> linear_index
   double  index_to_point (index_t ind) const {return xmin + ind * step;} 
   size_t  index_to_linear(index_t ind) const {return ind;}   
   //index_t point_to_index (double t) const {return (size_t)floor((t-x_min)/dt + 0.5);}    

   /// The wrapper for the mesh point
   struct mesh_point_t : arith_ops_by_cast<mesh_point_t, domain_pt_t  > {
    mesh_t const * m;  
    index_t index; 
    mesh_point_t( mesh_t const & mesh, index_t const & index_=0): m(&mesh), index(index_) {}
    void advance() { ++index;}
    operator domain_pt_t () const { return m->index_to_point(index);} 
   };

   /// Accessing a point of the mesh
   mesh_point_t operator[](index_t i) const { return mesh_point_t (*this,i);}

   /// Iterating on all the points...
   typedef  mesh_pt_generator<mesh_t> iterator;
   iterator begin() const { return iterator (this);}
   iterator end()   const { return iterator (this, true);}

   /// Mesh comparison
   bool operator == (mesh_t const & M) { return ((_dom == M._dom) && (L ==M.L) && (std::abs(x_min - M.x_min)<1.e-15) && (std::abs(x_max - M.x_max)<1.e-15));} 

   private:
   domain_t _dom;
   domain_pt_t xmin, xmax;
   size_t L; 
   double step;
  };
}}
#endif

