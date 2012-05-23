
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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

#ifndef TRIQS_ARRAYS_IMPL_COMMON_H
#define TRIQS_ARRAYS_IMPL_COMMON_H
#define TRIQS_ARRAYS_ALREADY_INCLUDED

// including python first remove some warning
#ifdef TRIQS_ARRAYS_WITH_PYTHON_SUPPORT
#include <boost/tuple/tuple_io.hpp>
#include "Python.h"
#endif

#include <exception> 

/// Maximum dimension of the arrays
#define ARRAY_NRANK_MAX 10

#ifdef __GNUC__
#define restrict __restrict__ 
#endif

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#define USE_STATIC_ASSERT
//#define USE_VARIADIC_TEMPLATES
#endif

#ifndef USE_STATIC_ASSERT
#include "boost/static_assert.hpp"
#define static_assert(X,MESS) BOOST_STATIC_ASSERT((X)) 
#endif

#include <assert.h>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

#define DISABLE_IF(Cond, Type) typename boost::disable_if< Cond, Type >::type
#define DISABLE_IFC(Cond, Type) typename boost::disable_if_c< Cond, Type >::type
#define ENABLE_IF(Cond, Type)  typename boost::enable_if< Cond, Type >::type
#define ENABLE_IFC(Cond, Type)  typename boost::enable_if_c< Cond, Type >::type

// Use Cblas
#define BOOST_NUMERIC_BINDINGS_BLAS_CBLAS 

namespace boost { namespace serialization { class access;}}

#include "../../utility/exceptions.hpp"
#include <sstream>
#define TRIQS_ARRAYS_THROW(s) { TRIQS_RUNTIME_ERROR<<s; } 
#define TRIQS_ARRAYS_CHECK_OR_THROW(Cond,Mess) {if (!(Cond)) {TRIQS_ARRAYS_THROW(Mess);}}

#ifdef TRIQS_ARRAYS_DEBUG
#define TRIQS_ARRAYS_DEBUG_CHECK(Cond,Error) TRIQS_ARRAYS_CHECK_OR_THROW(Cond,Error) 

#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#define TRIQS_ARRAYS_CHECK_IM_STORAGE_COMPAT
#define TRIQS_ARRAYS_ENFORCE_INIT_NAN_INF

#else 
#define TRIQS_ARRAYS_DEBUG_CHECK(Cond,Error) 
#endif

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_complex.hpp>

// Tags
namespace triqs { namespace arrays {

 namespace mpl=boost::mpl; //move to common

 struct use_default {};

 namespace Tag {	
  template<typename TAG, typename T> struct check: boost::is_base_of<TAG,T> {};
  struct indexmap_storage_pair{};
  struct expression{};
  struct expression_terminal{};
  struct array_algebra_expression_terminal{};
  struct matrix_algebra_expression_terminal{};
  struct vector_algebra_expression_terminal{};
  struct scalar_expression_terminal{};
  struct has_special_assign{}; 
  template <char C> struct has_special_infix {}; 
  struct has_immutable_array_interface{};
  struct no_init {}; struct nan_inf_init {}; struct default_init {};   
  
  struct h5_array_proxy {};

 }

 template <typename T> struct is_expression : Tag::check<Tag::expression,T> {};

 template<typename T> struct is_matrix_expr : Tag::check<Tag::matrix_algebra_expression_terminal,T> {}; 
 template<typename T> struct is_vector_expr : Tag::check<Tag::vector_algebra_expression_terminal,T> {}; 
 template<typename T> struct is_matrix_or_vector_expr : boost::mpl::or_<is_matrix_expr<T>, is_vector_expr<T> > {};

 template<typename T> struct has_immutable_array_interface : 
  boost::mpl::or_<
  Tag::check<Tag::has_immutable_array_interface,T>, 
  Tag::check<Tag::expression,T>, 
  Tag::check<Tag::indexmap_storage_pair,T> > {}; 

 template<typename T1,typename T2> struct has_immutable_array_interface2 :
  boost::mpl::and_<has_immutable_array_interface<T1>, has_immutable_array_interface<T2> > {}; 

 namespace Tag { struct array{}; struct array_view {}; struct C{}; struct Fortran{}; }
 template <typename T> struct is_array : Tag::check<Tag::array,T> {};
 template <typename T> struct is_array_view : Tag::check<Tag::array_view,T> {};
 template <typename T> struct is_array_or_view : boost::mpl::or_< is_array<T>, is_array_view<T> > {};

 namespace Tag { struct vector{}; struct vector_view {};}
 template <typename T> struct is_vector : Tag::check<Tag::vector,T> {};
 template <typename T> struct is_vector_view : Tag::check<Tag::vector_view,T> {};
 template <typename T> struct is_vector_or_view : boost::mpl::or_< is_vector<T>, is_vector_view<T> > {};

 namespace Tag { struct matrix_view {}; struct matrix {}; }
 template <typename T> struct is_matrix : Tag::check<Tag::matrix,T> {};
 template <typename T> struct is_matrix_view : Tag::check<Tag::matrix_view,T> {};
 template <typename T> struct is_matrix_or_view : boost::mpl::or_< is_matrix<T>, is_matrix_view<T> > {};

 template <typename T> struct has_special_assign : Tag::check<Tag::has_special_assign,T> {};

 // a lhs is either a immutableArray, or has a special assign and has a domain
 // It is understood that has_special_assign implies that has_a_domain
 // Do we want to add explicitely a has_domain ??
 template <typename T> struct is_array_assign_lhs : boost::mpl::or_< has_immutable_array_interface<T>, has_special_assign<T> >{};
 
 template <class T> struct is_value_class : boost::mpl::or_< is_array<T>, is_matrix<T>, is_vector<T> > {};
 template <class T> struct is_view_class : boost::mpl::or_< is_array_view<T>, is_matrix_view<T>, is_vector_view<T> > {};
 template <class T> struct is_value_or_view_class : boost::mpl::or_< is_value_class<T>, is_view_class<T> > {};

 // template<class S, class A> struct is_scalar_for : 
 //  boost::mpl::if_<boost::is_arithmetic<typename A::value_type > , boost::is_arithmetic<S>,  boost::is_same<S,typename A::value_type > > {};

 template <class S> struct is_scalar : boost::mpl::or_<boost::is_arithmetic<S > , boost::is_complex<S> > {};
 // too primitive ?
 template<class S, class A> struct is_scalar_for : 
  boost::mpl::if_<is_scalar<typename A::value_type > , is_scalar<S>,boost::is_same<S,typename A::value_type > >::type {};

 /// Makes a view
 template<typename A> typename A::view_type make_view(A const & x) { return typename A::view_type(x);}

 /// Makes a clone
 template<typename A> typename A::non_view_type make_clone(A const & x) { return typename A::non_view_type(x);}

 /// Is the data contiguous
 template<typename A> typename boost::enable_if<is_value_class<A>,bool>::type has_contiguous_data(A const &) {return true;}
 template<typename A> typename boost::enable_if<is_view_class<A>, bool>::type has_contiguous_data(A const & v){return v.indexmap().is_contiguous();}

 template< typename A> 
  typename boost::enable_if<is_view_class<A> >::type 
  resize_or_check_if_view ( A & a, typename A::shape_type const & sha) { 
   if (a.shape()!=sha) TRIQS_RUNTIME_ERROR<< "Size mismatch : view class shape = "<<a.shape() << " expected "<<sha;
  }

 template< typename A> 
  typename boost::enable_if<is_value_class<A> >::type 
  resize_or_check_if_view ( A & a, typename A::shape_type const & sha) { if (a.shape()!=sha) a.resize(sha); }

}}//namespace triqs::arrays
#endif

