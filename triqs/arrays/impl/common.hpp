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
#include <triqs/utility/first_include.hpp>
#include <triqs/clef.hpp>
#define TRIQS_ARRAYS_ALREADY_INCLUDED

/// Maximum dimension of the arrays
#define ARRAY_NRANK_MAX 10

#include <exception> 
#include <assert.h>
#include <triqs/utility/exceptions.hpp>
#include <sstream>

#include <boost/utility/enable_if.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/char.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <triqs/utility/compiler_details.hpp>
#include "./tags.hpp"
#include "./traits.hpp"
#include <triqs/utility/macros.hpp>

namespace boost { namespace serialization { class access;}}

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

namespace triqs { 

 typedef unsigned long long ull_t;
 
 /// Makes a view
 template<typename A> typename A::view_type make_view(A const & x) { return typename A::view_type(x);}

 /// Makes a clone
 template<typename A> typename A::non_view_type make_clone(A const & x) { return typename A::non_view_type(x);}

 namespace arrays {
  using triqs::make_clone;

  /// Is the data contiguous
  template<typename A> typename boost::disable_if<is_amv_value_or_view_class<A>,bool>::type has_contiguous_data(A const &) {return false;}
  template<typename A> typename boost::enable_if<is_amv_value_class<A>,bool>::type has_contiguous_data(A const &) {return true;}
  template<typename A> typename boost::enable_if<is_amv_view_class<A>, bool>::type has_contiguous_data(A const & v){return v.indexmap().is_contiguous();}

  template< typename A> 
   typename boost::enable_if<is_amv_view_class<A> >::type 
   resize_or_check_if_view ( A & a, typename A::shape_type const & sha) { 
    if (a.shape()!=sha) TRIQS_RUNTIME_ERROR<< "Size mismatch : view class shape = "<<a.shape() << " expected "<<sha;
   }

  template< typename A> 
   typename boost::enable_if<is_amv_value_class<A> >::type 
   resize_or_check_if_view ( A & a, typename A::shape_type const & sha) { if (a.shape()!=sha) a.resize(sha); }
 }}//namespace triqs::arrays
#endif

