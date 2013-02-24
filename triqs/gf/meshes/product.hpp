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
#ifndef TRIQS_GF_MESH_PRODUCT_H
#define TRIQS_GF_MESH_PRODUCT_H
#include "./mesh_tools.hpp"
#include "../domains/product.hpp"
namespace triqs { namespace gf {

 template<typename... Meshes> 
  struct mesh_product {
   typedef domain_product<typename Meshes::domain_t ... >  domain_t;
   typedef std::tuple<typename Meshes::index_t ... >       index_t; 
   typedef std::tuple<Meshes...>                           m_tuple_t;
   typedef typename domain_t::point_t                      domain_pt_t;

   static constexpr int dim = sizeof...(Meshes);

   mesh_product (Meshes const & ... meshes) : m_tuple(meshes...), _dom(meshes.domain()...) {}
   mesh_product () {}

   domain_t const & domain() const { return _dom;}

   size_t size() const {return _size(cint<dim-1>());}
   private: // implementation of size computation
   template<int N> size_t _size(cint<N> n) const { return _size(cint<N-1>()) * component(n).size();}
   size_t _size(cint<-1>) const { return 1;}
   public : 

   /// Each component of the mesh : mesh.component(one), mesh.component(zero)
   m_tuple_t const & components() const { return m_tuple;}
   
   template<int N> typename std::tuple_element<N,m_tuple_t>::type const & component(cint<N>) const { return std::get<N>(m_tuple);}
   template<int N> typename std::tuple_element<N,m_tuple_t>::type & component(cint<N>)             { return std::get<N>(m_tuple);}

   /// Conversions point <-> index <-> linear_index
   typename domain_t::point_t index_to_point(index_t const & ind) const { domain_pt_t res; index_to_point_impl(res, ind, cint<dim-1>()); return res; }
   private: // implementation
   template<int N> void index_to_point_impl( domain_pt_t & P, index_t const & ind, cint<N> n) const { 
    std::get<N>(P) = component(n).index_to_point(std::get<N>(ind));
    index_to_point_impl(P,ind,cint<N-1>());
   }
   void index_to_point_impl( domain_pt_t & P, index_t const & ind, cint<-1>) const {} 

   public:

   // index[0] + component[0].size * (index[1] + component[1].size* (index[2] + ....))
   size_t index_to_linear(index_t const & ind) const { size_t R=0; index_to_linear_impl(R,ind, cint<dim-1>()); return R; }
   private: // implementation
   template<int N> void index_to_linear_impl( size_t & R, index_t const & ind, cint<N> n) const { 
    R = component(n).index_to_linear(std::get<N>(ind)) + R * component(n).size();
    index_to_linear_impl(R,ind,cint<N-1>());
   }
   void index_to_linear_impl( size_t & R, index_t const & ind, cint<-1>) const {} 

   public:

   /// The wrapper for the mesh point
   struct mesh_point_t { 
    const mesh_product * m; 
    index_t index;

    /// The wrapper for one component of the mesh point
    template<int N, typename CastType> struct component : arith_ops_by_cast<component<N,CastType>, CastType>  { 
     mesh_point_t const * p; 
     component (mesh_point_t const & p_):p(&p_){}
     operator CastType() const { return std::get<N>(p->m->m_tuple).index_to_point(std::get<N>(p->index));}
    };
    template<int N> component<N,typename std::tuple_element<N,m_tuple_t>::type::domain_t::point_t> operator[](cint<N>) const { return *this;}

    mesh_point_t(mesh_product const & m_, index_t index_ ) : m(&m_), index(index_) {} 
    mesh_point_t(mesh_product const & m_)                  : m(&m_), index()    {} 

    typedef domain_pt_t cast_t;
    operator cast_t() const { return m->index_to_point(index);}

    // index[0] +=1; if index[0]==m.component[0].size() { index[0]=0; index[1] +=1; if  ....}  and so on until dim
    void advance() { advance_impl(cint<0>());}
    private: // implementation
    template<int N> void advance_impl(cint<N> n) { auto & i=std::get<N>(index); ++i; if (i==m->component(n).size()) {i=0;advance_impl(cint<N+1>());} }
    void advance_impl(cint<dim>) {}  
   };

   /// Accessing a point of the mesh
   mesh_point_t operator[](index_t i) const { return mesh_point_t (*this,i);}
   mesh_point_t operator()(typename Meshes::index_t ... i) const { return mesh_point_t (*this, std::make_tuple(i...));}

   /// Iterating on all the points...
   typedef  mesh_pt_generator<mesh_product> iterator;
   iterator begin() const { return iterator (this);}
   iterator end()   const { return iterator (this, true);}

   /// Mesh comparison
   friend bool operator == (mesh_product const & M1, mesh_product const & M2) { return M1.m_tuple==M2.m_tuple; }

   /// Write into HDF5
   friend void h5_write (h5::group fg, std::string subgroup_name, mesh_product const & m) {
    h5::group gr =  fg.create_group(subgroup_name);
    h5_write(gr,"domain",m.domain());
    m.h5_write_impl(gr,cint<0>());
   }

   private:
   template<int N> void h5_write_impl (h5::group gr, cint<N> n) {
     std::stringstream fs;fs <<"MeshComponent"<< N; 
     h5_write(gr,fs.str(), this->component(n));
     h5_write_impl(gr,cint<N+1>());
    }
   void h5_write_impl (h5::group gr, cint<dim>) {}

   /// Read from HDF5
   friend void h5_read  (h5::group fg, std::string subgroup_name, mesh_product & m){
    h5::group gr = fg.open_group(subgroup_name);
    h5_read(gr,"domain",m._dom);
    m.h5_read_impl(gr,cint<0>());
   }

   private:
   template<int N> void h5_read_impl (h5::group gr, cint<N> n) {
     std::stringstream fs; fs <<"MeshComponent"<< N; 
     h5_read(gr,fs.str(), this->component(n));
     h5_read_impl(gr,cint<N+1>());
    }
   void h5_read_impl (h5::group gr, cint<dim>) {}
 
   //  BOOST Serialization
   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
     ar & boost::serialization::make_nvp("domain",_dom);
    ser_impl(ar,version,cint<0>());
    }

   private:
   template<class Archive, int N> void ser_impl (Archive & ar, const unsigned int version, cint<N> n) {
     std::stringstream fs; fs <<"MeshComponent"<< N; 
     ar & boost::serialization::make_nvp(fs.str(),this->component(n));
     ser_impl(ar,version,cint<N+1>());
    }
   template<class Archive> void ser_impl (Archive & ar, const unsigned int version, cint<dim>) {}

   private:
   m_tuple_t  m_tuple;
   domain_t _dom;
  }; 

}}
#endif
