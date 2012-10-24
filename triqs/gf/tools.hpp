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
#ifndef TRIQS_GF_TOOLS_H
#define TRIQS_GF_TOOLS_H
#include <utility>
#include <boost/iterator/iterator_facade.hpp>
#include <triqs/clef/core.hpp>
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/utility/proto/tools.hpp>
#include <triqs/arrays/h5/simple_read_write.hpp>
#include <triqs/arrays/proto/matrix_algebra.hpp>
#include <triqs/arrays/proto/array_algebra.hpp>
#include "triqs/utility/complex_ops.hpp"
#include <triqs/utility/view_tools.hpp>

namespace triqs { namespace gf {

 namespace tqa= triqs::arrays;
 namespace mpl=boost::mpl;
 namespace proto=boost::proto;
 namespace bf = boost::fusion; 
 namespace tup = triqs::utility::proto;

 //------------------------------------------------------

 typedef std::complex<double> dcomplex; 

 enum statistic_enum {Boson,Fermion};

 struct freq_infty{}; // the point at infinity

 struct zero_t{}; struct one_t{}; struct two_t{}; struct three_t{};// to dispatch index_point in mesh...
 
 //------------------------------------------------------

 struct nothing {
  template<typename... Args> nothing(Args...) {} // takes anything, do nothing..
  nothing() {}
  typedef void has_view_type_tag;     // Idiom : ValueView  
  typedef nothing view_type;
  typedef nothing non_view_type;
   
  friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, nothing ) {}
  friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, nothing ) {}
 }; 

 //------------------------------------------------------

 // Derive from this object using CRTP to provide arithmetic operation by casting the final object to C 
 template<typename Derived, typename C> struct arith_ops_by_cast {};
#define IMPL_OP(OP)\
 template<typename D, typename C, typename T> \
 auto operator OP(arith_ops_by_cast<D,C> const & x, T y) -> decltype( std::declval<C>() OP y) {return C(static_cast<D const& >(x)) OP y;}\
 template<typename D, typename C, typename T> \
 auto operator OP(T y, arith_ops_by_cast<D,C> const & x) -> decltype (y OP std::declval<C>()) {return y OP C(static_cast<D const& >(x));}
 IMPL_OP(+); IMPL_OP(-); IMPL_OP(*); IMPL_OP(/);
#undef IMPL_OP

 //------------------------------------------------------

 template<typename MeshType>
  class mesh_pt_generator : 
   public boost::iterator_facade< mesh_pt_generator<MeshType>, typename MeshType::mesh_point_t const &, boost::forward_traversal_tag, 
   typename MeshType::mesh_point_t const & > {
    friend class boost::iterator_core_access;
    MeshType const * mesh;
    size_t u;
    typename MeshType::mesh_point_t pt;
    typename MeshType::mesh_point_t const & dereference() const { return pt;}
    bool equal(mesh_pt_generator const & other) const { return ((mesh == other.mesh) && (other.u==u) );}
    public:
    mesh_pt_generator( MeshType const * m=NULL, bool atEnd = false): mesh(m), u(atEnd ? m->size(): 0), pt(*m) {}
    //~mesh_pt_generator() { std::cout  << "deleting generator" << std::endl; }
    void increment() { ++u; pt.advance(); }
    bool at_end() const { return (u>=mesh->size()-1);}
    typename MeshType::domain_t::point_t to_point() const { return pt;}    
   };
 //------------------------------------------------------

 /// The wrapper for the mesh point
 template<typename MeshType>
  struct mesh_point_d1 : arith_ops_by_cast<mesh_point_d1<MeshType>, typename MeshType::domain_t::point_t  > {
   typedef MeshType mesh_t;
   mesh_t const * m;  
   typename mesh_t::index_t index; 
   mesh_point_d1( mesh_t const & mesh, typename mesh_t::index_t const & index_=0): m(&mesh), index(index_) {}
   void advance() { ++index;}
   operator typename MeshType::domain_t::point_t () const { return m->index_to_point(index);} 
  };

 //------------------------------------------------------

 template<typename Domain1, typename Domain2> 
  struct domain_product_2 { 
   typedef std::tuple<typename Domain1::point_t, typename Domain2::point_t> point_t;
   Domain1 d1; 
   Domain2 d2;
   domain_product_2(Domain1 const & dd1, Domain2 const & dd2): d1(dd1), d2(dd2) {}
  };

 //------------------------------------------------------

 template<typename Mesh0, typename Mesh1> 
  struct mesh_product_2 {
   typedef domain_product_2<typename Mesh0::domain_t, typename Mesh1::domain_t>  domain_t;
   typedef std::tuple<typename Mesh0::index_t, typename Mesh1::index_t>          index_t; 

   std::tuple<Mesh0,Mesh1> m;
   mesh_product_2 (Mesh0 const & m0, Mesh1 const & m1) : m(m0,m1), _dom(domain_t(m0.domain(),m1.domain())){}
   mesh_product_2 () {}

   domain_t const & domain() const { return _dom;}

   size_t size() const {return std::get<0>(m).size() * std::get<1>(m).size();}

   /// Conversions point <-> index <-> linear_index
   typename domain_t::point_t index_to_point(index_t const & ind) const { 
    return std::make_tuple( std::get<0>(m).index_to_point(std::get<0>(ind)), std::get<1>(m).index_to_point(std::get<1>(ind)));
   }

   size_t index_to_linear(index_t const & ind) const { 
    return std::get<0>(m).index_to_linear(std::get<0>(ind)) + std::get<1>(m).index_to_linear(std::get<1>(ind))*std::get<0>(m).size();
   }

   /// The wrapper for the mesh point
   struct mesh_point_t { 
    const mesh_product_2 * m; index_t index;

    /// The wrapper for one component of the mesh point
    template<int N, typename CastType> struct component : arith_ops_by_cast<component<N,CastType>, CastType>  { 
     mesh_point_t const * p; 
     component (mesh_point_t const & p_):p(&p_){}
     operator CastType() const { return std::get<N>(p->m->m).index_to_point(std::get<N>(p->index));}
    };
    component<0,typename Mesh0::domain_t::point_t> _0() const { return *this;}
    component<1,typename Mesh1::domain_t::point_t> _1() const { return *this;}

    mesh_point_t(mesh_product_2 const & m_, index_t index_ ) : m(&m_), index(index_) {} 
    mesh_point_t(mesh_product_2 const & m_)                  : m(&m_), index()    {} 

    void advance() { ++std::get<0>(index); if (std::get<0>(index)==std::get<0>(m->m).size()) { std::get<0>(index)=0; ++std::get<1>(index);}}
    typedef typename domain_t::point_t cast_t;
    operator cast_t() const { return m->index_to_point(index);}
   };

   /// Accessing a point of the mesh
   mesh_point_t operator[](index_t i) const { return mesh_point_t (*this,i);}
   mesh_point_t operator()(typename Mesh0::index_t i, typename Mesh1::index_t j) const { return mesh_point_t (*this, std::make_tuple(i,j));}

   /// Iterating on all the points...
   typedef  mesh_pt_generator<mesh_product_2> iterator;
   iterator begin() const { return iterator (this);}
   iterator end()   const { return iterator (this, true);}

   /// Mesh comparison
   friend bool operator == (mesh_product_2 const & M1, mesh_product_2 const & M2) { return ( (std::get<0>(M1.m)==std::get<0>(M2.m)) && (std::get<1>(M1.m)==std::get<1>(M2.m))); }

   private:
   domain_t _dom;
  }; 

}}
#endif
