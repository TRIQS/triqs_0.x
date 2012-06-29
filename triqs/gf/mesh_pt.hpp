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
#ifndef TRIQS_GF_MESH_PT_H
#define TRIQS_GF_MESH_PT_H
#include "domains.hpp"

namespace triqs { namespace gf { namespace meshes { 

 namespace details { 
  // a trick to avoid putting T() in the typeof type deduction ! This code is NEVER used 
  template<class T> typename boost::unwrap_reference<T>::type pdc() { 
   typename boost::unwrap_reference<T>::type * x= NULL; assert(0); return *x; 
  }

  // a tool to compute the return type of basic ops...
  template<typename A, typename B> struct _add_ { typedef BOOST_TYPEOF_TPL( pdc <A>() + pdc<B>()) type;};
  template<typename A, typename B> struct _sub_ { typedef BOOST_TYPEOF_TPL( pdc <A>() - pdc<B>()) type;};
  template<typename A, typename B> struct _mul_ { typedef BOOST_TYPEOF_TPL( pdc <A>() * pdc<B>()) type;};
  template<typename A, typename B> struct _div_ { typedef BOOST_TYPEOF_TPL( pdc <A>() / pdc<B>()) type;};
 }

 template<typename MeshType> struct mesh_pt {
  MeshType const & m;
  typename MeshType::index_type i;
  mesh_pt( MeshType const & mesh, typename MeshType::index_type const & index): m(mesh), i(index) {}
  typedef typename MeshType::domain_type::embedded_point_type cast_type;
  operator cast_type() const { return m.embed(i);} // cast into the element type of the domain (e.g. real time, real frequency).
  cast_type casted() const { return m.embed(i);}
  // I need to explicitely redefine the 4 basic operations...
#define IMPL_OP(TRA, OP)\
  template<typename T> friend typename details::TRA<cast_type,T>::type operator OP ( mesh_pt const & x, T y) { return cast_type(x) OP y;}\
  template<typename T> friend typename details::TRA<T,cast_type>::type operator OP (T y, mesh_pt const & x)  { return y OP cast_type(x);}
  IMPL_OP(_add_,+);
  IMPL_OP(_sub_,-);
  IMPL_OP(_mul_,*);
  IMPL_OP(_div_,/);
#undef IMPL_OP
 };

 template<typename MeshType> 
  mesh_pt<MeshType> make_mesh_pt(MeshType const & m, typename MeshType::index_type const & i){ return mesh_pt<MeshType>(m,i);}

}}}

#endif




