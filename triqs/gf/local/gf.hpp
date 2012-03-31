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

#ifndef TRIQS_GF_LOCAL_GF_H
#define TRIQS_GF_LOCAL_GF_H 

#include <triqs/lazy/core.hpp>
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>

// general
namespace triqs { namespace gf { 
 
 namespace tqa= triqs::arrays;
 namespace tql= triqs::lazy;
  
 template<typename LHS, typename RHS, typename Enable=void > struct assignment_impl;
}}


namespace triqs { namespace gf { namespace local {

//template <typename T> struct view_type_of { typedef typename T::view_type type;};
//template <> struct view_type_of<void> { typedef void type;};

struct nothing { 

 typedef nothing view_type;
 template<typename A, typename B> view_type slice(A a, B b) const { return nothing();}

};

using tqa::range;

 template<typename DomainType, bool IsView, typename AuxType = nothing> class gf {

 public:

  typedef double value_element_type;
  typedef DomainType domain_type;
  typedef tqa::matrix_view<element_value_type, Option::Fortran>       result_type;
  typedef tqa::matrix_view<const element_value_type, Option::Fortran> const_result_type;

  typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
  typedef gf<DomainType,true>   view_type;
  typedef gf<DomainType,false>  non_view_type;

  typedef std::vector<string> indices_type;

  //typedef AuxType aux_type;
  //typedef typename view_type_of<aux_type>::type aux_view_type; // correct this
  //typedef mpl::if_< boost::is_same<aux_type,void>, void, typename aux_type::typename view_type_of<aux_type>::type aux_view_type; // correct this
  //typedef typename mpl::if_<IsView, aux_view_type, aux_type>::type  aux_type;
  typedef typename AuxType::view_type                              aux_view_type; 
  typedef typename mpl::if_<IsView, aux_view_type, AuxType>::type  aux_type;

  typedef tqa::array      <value_element_type,3>                               data_non_view_type;
  typedef tqa::array_view <value_element_type,3>                               data_view_type;
  typedef typename mpl::if_<IsView, data_view_type, data_non_view_type>::type  data_type;
 
  typedef typename DomainType::index_type arg0_type;

 protected:
  mesh_type mesh;
  data_type data;
  aux_type aux;
  indices_type indices_left, indices_right;

 public:

  // not valid for the view, need a disabler
  template<typename MeshConstructorArgType> 
   gf (size_t N1, size_t N2, MeshConstructorArgType const & mesh_, indices_type const & indices_left_, indices_type const & indices_right_, 
     typename boost::disable_if_c<IsView>::type * dummy =NULL ) : 
    mesh(mesh_), data(N1,N2,mesh_.len()), indices_left(indices_left_), indices_right(indices_right_) {
     // init aux
    }

  gf (mesh_type const & mesh_, data_view_type const & data_, aux_view_type const & aux_, indices_type const & indices_left_, indices_type const & indices_right_) : 
   data(data_), mesh(mesh_), aux(aux_), indices_left(indices_left_), indices_right(indices_right_) {}

  gf(gf const & x): mesh(x.mesh), data(x.data), aux(x.aux), indices_left(x.indices_left), indices_right(x.indices_right) {}

  // view from regular, and vice versa
  template<typename GfType> gf(GfType const & x): mesh(x.mesh), data(x.data), aux(x.aux), indices_left(x.indices_left), indices_right(x.indices_right) {}

  //template<typename T> typename boost::disable_if< tql::is_lazy<T>, result_type>::type operator() (T const & x) { return data(range(),range(),x);}
  //template<typename T> typename boost::disable_if< tql::is_lazy<T>, const_result_type>::type operator() (T const & x) const { return data(range(),range(),x);}

  triqs::arrays::mini_vector<size_t,2> result_shape() const { return data.shape();}
  
  result_type       operator()       ( arg0_type const & x) { return data(range(),range(),x);}
  const_result_type operator() const ( arg0_type const & x) { return data(range(),range(),x);}

  aux_view_type       operator()       ( domain::infty const & x) { return aux;}
  const aux_view_type operator() const ( domain::infty const & x) { return aux;}

  TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,gf_view);
  //TRIQS_LAZY_ADD_LAZY_CALL_WITH_COPY(1,gf);

  template<typename A, typename B> view_type slice(A a, B b) { 
   return view_type( mesh, data (range(a), range(b), range()), aux.slice(a,b), _myslice(indices_left,range(a)), _myslice(indices_right,range(b)) ); 
  }
  template<typename A, typename B> const view_type slice(A a, B b) const { 
   return view_type( mesh, data (range(a), range(b), range()), aux.slice(a,b), _myslice(indices_left,range(a)), _myslice(indices_right,range(b)) ); 
  }

  private:
  //template <typename D, bool IV, typename A, typename T1, typename T2> gf<D,true,A> _myslice( gf<D,IV,A> const &X, T1 a, T2 b) { return X.slice(a,b);}
  //template <typename X, typename T1, typename T2> typename X::view_type _myslice( X const &x, T1 a, T2 b) { return x.slice(a,b);}
  //template <typename T1, typename T2> void _myslice( void, T1 a, T2 b) {}
  std::vector<string> _myslice(std::vector<string> const & V, range const & R) { }

  public:

  //view_type view()             { return view_type(mesh,data,aux,indices_left, indices_right);}
  //const view_type view() const { return view_type(mesh,data,aux,indices_left, indices_right);}

  template<typename RHS> // specialize for various RHS ( fourier_impl, other gf, etc....)
   gf & operator = (RHS const & rhs) { triqs::gf::assignment_impl<gf,RHS>::invoke(*this,rhs); return *this; } 

  // lazy_assignable
  template<typename F>
   void set_from_function(F f) { const size_t Nmax = data.shape()[2]; for (size_t u=0; u<Nmax; ++u) data(range(),range(),u) = F(u);}

  /// Save the Green function in i omega_n (as 2 columns).
  void save(string file,  bool accumulate=false) const;

  /// Load the GF
  void load(string file);

  /// HDF5 saving ....

 };

 // gf.impl.hpp

 // implementation of expression template as function of omega...
 // keep the domain

 // 
 // integrate + enable ??
 template< typename GfType> struct integrate_impl; 
 template< typename GfType> typename GfType::result_type integrate(GfType const & x) { return integrate_impl<GfType>::invoke(x);}
 // + lazy overload + work on the concept only !!
 // 
 // density + enable ??
 template< typename GfType> struct density_impl; 
 template< typename GfType> typename GfType::result_type density(GfType const & x) { return density_impl<GfType>::invoke(x);}
 //+ lazy overload.

 // 1 fichier de specialization par type de GF
 //

 // rewrite the immutable_concept : call, slice, domain

 // implementation of save/load

}}

namespace lazy { 

 // overload of the free functions



}

}
#endif




