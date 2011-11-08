
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

#ifndef TRIQS_ARRAYS_EXPRESSION_MINMAX_H
#define TRIQS_ARRAYS_EXPRESSION_MINMAX_H
#include <algorithm>
#include <functional>
#include "./fold.hpp"

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000  + __GNUC_MINOR__ * 100  + __GNUC_PATCHLEVEL__)
#endif

namespace triqs { namespace arrays {

// fix for gcc <4.4. Remove it to go for C++11
// clang defines __GNUC__ , so we explicitely remove this case.

#if (!defined(__clang__)) and defined (GCC_VERSION) and (GCC_VERSION < 40400)

 template <class T> T mymax ( T a, T b ) { return std::max(a,b); }
 template <class T> T mymin ( T a, T b ) { return std::min(a,b); }

 template<class A>
  typename A::value_type max_element(A const &a) { 
   return fold ( mymax<typename A::value_type  >)  ( a, a[typename A::domain_type::index_value_type()]);// init = first element
  }

 template<class A>
  typename A::value_type min_element(A const &a) { 
   return fold ( mymin<typename A::value_type  >)  ( a, a[typename A::domain_type::index_value_type()]);
  }

#else

 template<class A>
  typename A::value_type max_element(A const &a) { 
   return fold ( std::max<typename A::value_type const & >)  ( a, a[typename A::domain_type::index_value_type()]);// init = first element
  }

 template<class A>
  typename A::value_type min_element(A const &a) { 
   return fold ( std::min<typename A::value_type const & >)  ( a, a[typename A::domain_type::index_value_type()]);
  }

#endif

}}//namespace triqs::arrays 

#endif


