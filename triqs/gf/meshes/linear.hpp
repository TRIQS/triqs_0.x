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
#include "./mesh_tools.hpp"
namespace triqs { namespace gf { 

 // Three possible meshes
 enum mesh_kind { half_bins, full_bins, without_last };

 template<typename Domain>
  struct linear_mesh {

   typedef Domain domain_t;
   typedef size_t index_t; 
   typedef typename domain_t::point_t  domain_pt_t;

   linear_mesh (domain_t const & dom, double a, double b, size_t n_pts, mesh_kind mk) :
     _dom(dom), L(n_pts), a_pt(a), b_pt(b), meshk(mk) {
     switch(mk) {
       case half_bins: del = (b-a)/L; xmin = a+0.5*del; break;
       case full_bins: del = (b-a)/(L-1); xmin = a; break;
       case without_last: del = (b-a)/L; xmin = a; break;
     }
     xmax = xmin + del*(L-1);
   }

   linear_mesh (domain_t && dom, double a, double b, size_t n_pts, mesh_kind mk) :
     _dom(dom), L(n_pts), a_pt(a), b_pt(b), meshk(mk) {
     switch(mk) {
       case half_bins: del = (b-a)/L; xmin = a+0.5*del; break;
       case full_bins: del = (b-a)/(L-1); xmin = a; break;
       case without_last: del = (b-a)/L; xmin = a; break;
     }
     xmax = xmin + del*(L-1);
   }

   linear_mesh () : _dom(), L(0), a_pt(0), b_pt(0), xmin(0), xmax(0), del(0), meshk(half_bins) {}

   domain_t const & domain() const { return _dom;}
   size_t size() const { return L; }
   double delta() const { return del; }
   double x_max() const { return xmax; }
   double x_min() const { return xmin; }
   mesh_kind kind() const { return meshk; }

   /// Conversions point <-> index <-> linear_index
   domain_pt_t  index_to_point (index_t ind) const {return embed(xmin + ind * del, mpl::bool_<boost::is_base_of<std::complex<double>, domain_pt_t>::value >()) ;}
   private : // multiply by I is the type is a complex ....
   domain_pt_t embed( double x, mpl::bool_<false> ) const { return x;}
   domain_pt_t embed( double x, mpl::bool_<true> ) const { return std::complex<double>(0,x);}
   public : 

   size_t  index_to_linear(index_t ind) const {return ind;}   
 
   /// The wrapper for the mesh point
   class mesh_point_t : tag::mesh_point, public arith_ops_by_cast<mesh_point_t, domain_pt_t  > {
    linear_mesh const * m;  
    index_t _index; 
    public:
    mesh_point_t( linear_mesh const & mesh, index_t const & index_): m(&mesh), _index(index_) {}
    void advance() { ++_index;}
    operator domain_pt_t () const { return m->index_to_point(_index);} 
    size_t linear_index() const { return _index;}
    size_t index() const { return _index;}
    bool at_end() const { return (_index == m->size());}
    void reset() {_index =0;}
   };

   /// Accessing a point of the mesh
   mesh_point_t operator[](index_t i) const { return mesh_point_t (*this,i);}

   /// Iterating on all the points...
   typedef  mesh_pt_generator<linear_mesh> iterator;
   iterator begin() const { return iterator (this);}
   iterator end()   const { return iterator (this, true);}

   /// Mesh comparison
   bool operator == (linear_mesh const & M) const { return ((_dom == M._dom) && (size() ==M.size()) && (std::abs(xmin - M.xmin)<1.e-15) && (std::abs(xmax - M.xmax)<1.e-15));} 

   /// Write into HDF5
   friend void h5_write (h5::group fg, std::string subgroup_name, linear_mesh const & m) {
    h5::group gr =  fg.create_group(subgroup_name);
    int k;
    switch(m.meshk) {
       case half_bins: k=0; break;
       case full_bins: k=1; break;
       case without_last: k=2; break;
    }
    h5_write(gr,"domain",m.domain());
    h5_write(gr,"min",m.a_pt);
    h5_write(gr,"max",m.b_pt);
    h5_write(gr,"size",m.size());
    h5_write(gr,"kind",k);
   }

   /// Read from HDF5
   friend void h5_read  (h5::group fg, std::string subgroup_name, linear_mesh & m){
    h5::group gr = fg.open_group(subgroup_name);
    typename linear_mesh::domain_t dom;
    double a,b;
    size_t L; 
    int k;
    mesh_kind mk;
    h5_read(gr,"domain",dom);
    h5_read(gr,"min",a);
    h5_read(gr,"max",b);
    h5_read(gr,"size",L);
    h5_read(gr,"kind",k);
    switch(k) {
       case 0: mk = half_bins; break;
       case 1: mk = full_bins; break;
       case 2: mk = without_last; break;
    }
    m = linear_mesh(std::move(dom), a, b, L, mk);
   }

   //  BOOST Serialization
   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
     ar & boost::serialization::make_nvp("domain",_dom);
     ar & boost::serialization::make_nvp("a_pt",a_pt);
     ar & boost::serialization::make_nvp("b_pt",b_pt);
     ar & boost::serialization::make_nvp("xmin",xmin);
     ar & boost::serialization::make_nvp("xmax",xmax);
     ar & boost::serialization::make_nvp("del",del);
     ar & boost::serialization::make_nvp("size",L);
     ar & boost::serialization::make_nvp("kind",meshk);
    }

   private:
   domain_t _dom;
   size_t L; 
   double a_pt, b_pt;
   double xmin, xmax;
   double del;
   mesh_kind meshk;
  };

 /// Approximation of a point of the domain by a mesh point  
 template<typename D>
  std::tuple<bool, size_t, double>  windowing ( linear_mesh<D> const & mesh, typename D::point_t const & x) { 
   double a = (x - mesh.x_min())/mesh.delta();
   long i = floor(a);
   bool in = (! ((i<0) || (i>=long(mesh.size())-1)));
   double w = a-i;
   //   std::cerr  << " window "<< i << " "<< in << "  "<< w<< std::endl ;
   return std::make_tuple(in, (in ? size_t(i) : 0),w);
 }

}}
#endif

