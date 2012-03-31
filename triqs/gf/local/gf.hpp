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

namespace triqs { namespace gf { namespace local {

 namespace tqa= triqs::arrays;
 namespace tql= triqs::lazy;
 using tqa::range;

 template<typename DomainType, bool IsView> class gf;

 // implementation details declaration
 namespace impl {
  template<typename LHS, typename RHS, typename Enable=void > struct assignment;
  inline std::vector<string> slice_vector(std::vector<string> const & V, range const & R);
  struct nothing; 
  template<typename DomainType> struct tail_type_from_domain;
  template<typename T, typename F> struct assign_to_F_of_infty;
 }

 namespace impl { // separate ?
  struct nothing { 
   typedef nothing view_type;
   template<typename A, typename B> view_type slice(A a, B b) const { return nothing();}
   nothing(...){}
   void operator=(...) {}
  };

  template<typename DomainType> struct tail_type_from_domain { typedef nothing type;}
  template<> struct tail_type_from_domain<domains::matsubara_freq> { typedef gf<domains::tail,false> type;}

  template<typename T, typename F> struct assign_to_F_of_infty { static void invoke(T & lhs, F const & f) { lhs = F(domains::infty);} };
  template<typename F> struct assign_to_F_of_infty<nothing,F> { static void invoke(T & lhs, F const & f) {} };
 }

 // the class 
 template<typename DomainType, bool IsView> class gf {

  public:

   typedef double value_element_type;
   typedef DomainType domain_type;
   typedef std::vector<string> indices_type;

   typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
   typedef gf<DomainType,true>   view_type;
   typedef gf<DomainType,false>  non_view_type;

   typedef typename tail_type_from_domain<DomainType>::type                      tail_non_view_type; 
   typedef typename tail_non_view_type::view_type                                tail_view_type; 
   typedef typename mpl::if_<IsView, tail_view_type, tail_non_view_type>::type   tail_type;

   typedef tqa::array      <value_element_type,3>                               data_non_view_type;
   typedef tqa::array_view <value_element_type,3>                               data_view_type;
   typedef typename mpl::if_<IsView, data_view_type, data_non_view_type>::type  data_type;

   typedef typename DomainType::index_type arg0_type;
   typedef tqa::matrix_view<element_value_type, Option::Fortran>       result_type;
   typedef tqa::matrix_view<const element_value_type, Option::Fortran> const_result_type;

  protected:
   mesh_type mesh;
   data_type data;
   tail_type tail;
  public:

   const indices_type indices_left, indices_right;

   mesh_type const & mesh() const { return mesh;}

   data_view_type data_view() { return data;}
   const data_view_type data_view() const { return data;}

   tail_view_type tail_view() { return tail;}
   const tail_view_type tail_view() const { return tail;}

   template<typename S> // not valid for the view, need a template and a disabler 
    gf (S N1, S N2, mesh_type const & mesh_, indices_type const & indices_left_, indices_type const & indices_right_, 
      typename boost::disable_if_c<IsView>::type * dummy =NULL ) : 
     mesh(mesh_), data(N1,N2,mesh_.len()),tail(mesh_.mesh_tail),  indices_left(indices_left_), indices_right(indices_right_) {
     }

   /*  // not valid for the view, need a disabler
       template<typename S, typename Enable = typename boost::disable_if_c<IsView>::type > 
       gf (S N1, S N2, mesh_type const & mesh_, indices_type const & indices_left_, indices_type const & indices_right_):
       mesh(mesh_), data(N1,N2,mesh_.len()), tail(mesh_.mesh_tail), indices_left(indices_left_), indices_right(indices_right_) {
       }
       */

   gf (mesh_type const & mesh_, data_view_type const & data_, tail_view_type const & tail_, indices_type const & indices_left_, indices_type const & indices_right_) : 
    data(data_), mesh(mesh_), tail(tail_), indices_left(indices_left_), indices_right(indices_right_) {}

   gf(gf const & x): mesh(x.mesh), data(x.data), tail(x.tail), indices_left(x.indices_left), indices_right(x.indices_right) {}

   template<typename GfType> gf(GfType const & x): mesh(x.mesh), data(x.data), tail(x.tail), indices_left(x.indices_left), indices_right(x.indices_right) {}

   triqs::arrays::mini_vector<size_t,2> result_shape() const { return data.shape();}

   result_type       operator()       ( arg0_type const & x) { return data(range(),range(),x);}
   const_result_type operator() const ( arg0_type const & x) { return data(range(),range(),x);}

   tail_view_type       operator()       ( domain::infty const & x) { return tail;}
   const tail_view_type operator() const ( domain::infty const & x) { return tail;}

   //TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);
   TRIQS_LAZY_ADD_LAZY_CALL_WITH_COPY(1,view_type);

   template<typename A, typename B> view_type slice(A a, B b) { 
    return view_type( mesh, data (range(a), range(b), range()), tail.slice(a,b), impl::slice_vector(indices_left,range(a)), impl::slice_vector(indices_right,range(b)) ); 
   }
   template<typename A, typename B> const view_type slice(A a, B b) const { 
    return view_type( mesh, data (range(a), range(b), range()), tail.slice(a,b), impl::slice_vector(indices_left,range(a)), impl::slice_vector(indices_right,range(b)) ); 
   }

   template<typename RHS> // specialize for various RHS ( fourier_impl, other gf, etc....)
    gf & operator = (RHS const & rhs) { triqs::gf::impl::assignment<gf,RHS>::invoke(*this,rhs); return *this; } 

   // lazy_assignable
   template<typename F>
    void set_from_function(F f) { 
     const size_t Nmax = data.shape()[2]; for (size_t u=0; u<Nmax; ++u) data(range(),range(),u) = F(u);
     impl::assign_to_F_of_infty<tail_type,F>::invoke(tail,f);
    }

   /// Save the Green function in i omega_n (as 2 columns).
   void save(string file,  bool accumulate=false) const {}

   /// Load the GF
   void load(string file){}

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
 // rewrite the immutable_concept : call, slice, domain


}}}
#endif




