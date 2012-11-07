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
#ifndef TRIQS_UTILITY_VIEWTOOLS_H
#define TRIQS_UTILITY_VIEWTOOLS_H
#include <functional>

namespace triqs { 

 template<typename T, typename Void =void> struct view_type_if_exists_else_type {typedef T type;}; 
 template<typename T> struct view_type_if_exists_else_type<T, typename T::has_view_type_tag> {typedef typename T::view_type type;}; 
 
 template<typename T, typename Void =void> struct const_view_type_if_exists_else_type {typedef T type;}; 
 template<typename T> struct const_view_type_if_exists_else_type<T, typename T::has_view_type_tag> {typedef const typename T::view_type type;}; 
 // template<typename T> struct const_view_type_if_exists_else_type<const T, typename T::has_view_type_tag> {typedef const typename T::view_type type;}; 

 // replacement of std::plus for views ...
 template <class T> struct add_views : std::binary_function <T,T,T> {
  T operator() (const T& x, const T& y) const
  { typename T::non_view_type r(x); r =r + y; return r;}
 };

 // This is a little proxy, the only class allowed to build a view empty !
 template<typename V> class view_proxy {
  V _v;
  public:
  view_proxy() {}
  view_proxy(V const &v) :_v(v){}
  void rebind( V const & v) { _v.rebind(v);}
  V & operator() () { return _v;}
  V const & operator() () const { return _v;}
  template<typename V2> void operator <<( V2 && v) { _v = v;} 
 };

}//namespace triqs
#endif

