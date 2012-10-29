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

   linear_mesh (domain_t && dom, double x_min, double x_max, size_t ns) : _dom(dom),L(ns), xmin(x_min), xmax(x_max), _step( (x_max - x_min) /(ns-1)) {}
   linear_mesh () : _dom(),L(0){}

   domain_t const & domain() const { return _dom;}
   size_t size() const {return L;}
   double delta() const { return _step;}
   double x_max() const { return xmax;}
   double x_min() const { return xmin;}

   /// Conversions point <-> index <-> linear_index
   domain_pt_t  index_to_point (index_t ind) const {return embed( xmin + ind * _step, mpl::bool_<boost::is_base_of<std::complex<double>, domain_pt_t>::value >()) ;}
   private : // multiply by I is the type is a complex ....
   domain_pt_t embed( double x, mpl::bool_<false> ) const { return x;}
   domain_pt_t embed( double x, mpl::bool_<true> ) const { return std::complex<double>(0,x);}
   public : 

   size_t  index_to_linear(index_t ind) const {return ind;}   

   std::tuple<index_t,double,bool> closest_point (double t) const {
    size_t i = floor((t-xmin)/_step + 0.5); 
    double delta = (t-xmin) - _step * i;
    return std::make_tuple(i,delta, ((i>=0) && (i< L)));
   }

   /// The wrapper for the mesh point
   struct mesh_point_t : arith_ops_by_cast<mesh_point_t, domain_pt_t  > {
    linear_mesh const * m;  
    index_t index; 
    mesh_point_t( linear_mesh const & mesh, index_t const & index_=0): m(&mesh), index(index_) {}
    void advance() { ++index;}
    operator domain_pt_t () const { return m->index_to_point(index);} 
   };

   /// Accessing a point of the mesh
   mesh_point_t operator[](index_t i) const { return mesh_point_t (*this,i);}

   /// Iterating on all the points...
   typedef  mesh_pt_generator<linear_mesh> iterator;
   iterator begin() const { return iterator (this);}
   iterator end()   const { return iterator (this, true);}

   /// Mesh comparison
   bool operator == (linear_mesh const & M) const { return ((_dom == M._dom) && (L ==M.L) && (std::abs(xmin - M.xmin)<1.e-15) && (std::abs(xmax - M.xmax)<1.e-15));} 

   /// Write into HDF5
   friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, linear_mesh const & m) {
    tqa::h5::group_or_file gr =  fg.create_group(subgroup_name);
    h5_write(gr,"Domain",m.domain());
    // COMPLEX NOT H5 READ/WRITE !! FORGOTTEN ...
    h5_write(gr,"min",m.xmin);
    h5_write(gr,"max",m.xmax);
    h5_write(gr,"size",m.size());
   }

   /// Read from HDF5
   friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, linear_mesh & m){
    tqa::h5::group_or_file gr = fg.open_group(subgroup_name);
    typename linear_mesh::domain_t dom;
    typename linear_mesh::domain_pt_t xmin,xmax;
    size_t Nmax; 
    h5_read(gr,"Domain",m._dom);
    h5_read(gr,"min",xmin);
    h5_read(gr,"max",xmax);
    h5_read(gr,"size",Nmax);
    m = linear_mesh(std::move(dom), xmin,xmax,Nmax);
   }

   private:
   domain_t _dom;
   double xmin, xmax;
   size_t L; 
   double _step;
  };
}}
#endif

