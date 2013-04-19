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
#ifndef TRIQS_ARRAYS_ALGORITHMS_H
#define TRIQS_ARRAYS_ALGORITHMS_H
#include<triqs/utility/first_include.hpp>
#include <algorithm>
#include <functional>
#include "./functional/fold.hpp"
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

#define TRIQS_FIX_CLANG30_linux

namespace triqs { namespace arrays {

 // get_first_element(a) returns .. the first element, i.e. a(0,0,0...)
#define ZERO(z, n, text) 0
#define IMPL(z, N, unused)\
 template<typename ArrayType> typename std::enable_if<ArrayType::domain_type::rank==N,typename ArrayType::value_type>::type get_first_element(ArrayType const & A) \
 { return A(BOOST_PP_ENUM(N,ZERO,0));}
 BOOST_PP_REPEAT_FROM_TO(1,ARRAY_NRANK_MAX , IMPL, nil);
#undef IMPL
#undef ZERO

 template<class A>
  typename A::value_type max_element(A const &a) { 
#ifdef TRIQS_FIX_CLANG30_linux
   typedef typename A::value_type const & T;
   return fold ( static_cast<T(*)(T,T)>(std::max<T>))  ( a, get_first_element(a)); //a[typename A::domain_type::index_value_type()]);
#else  
   return fold ( std::max<typename A::value_type const & >)  ( a, get_first_element(a)); //a[typename A::domain_type::index_value_type()]);
#endif
  }

 template<class A>
  typename A::value_type min_element(A const &a) { 
#ifdef TRIQS_FIX_CLANG30_linux
 typedef typename A::value_type const & T;
   return fold ( static_cast<T(*)(T,T)>(std::min<T>))  ( a, get_first_element(a)); //a[typename A::domain_type::index_value_type()]);
#else 
   return fold ( std::min<typename A::value_type const & >)  ( a,  get_first_element(a)); //a[typename A::domain_type::index_value_type()]);
#endif
  }

 template <class A>
  typename A::value_type sum(A const & a) { return fold ( std::plus<typename A::value_type>())  (a); }

 template <class A>
  typename A::value_type prod(A const & a) { return fold ( std::multiplies<typename A::value_type>())  (a,1); }


}}//namespace triqs::arrays 

#endif


