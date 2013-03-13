/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by O. Parcollet
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
#ifndef TRIQS_APPLY_ON_TUPLE_H
#define TRIQS_APPLY_ON_TUPLE_H
#include<triqs/utility/macros.hpp>
#include <tuple>
namespace triqs {

 template<int pos,  typename T> struct apply_on_tuple_impl {
  template<typename F, typename ... Args>
   auto operator()(F && f, T const & t, Args && ... args)
   DECL_AND_RETURN( apply_on_tuple_impl<pos-1,T>()(std::forward<F>(f),t, std::get<pos>(t), std::forward<Args>(args)...));
 };

 template<typename T> struct apply_on_tuple_impl<-1,T> {
  template<typename F, typename ... Args>
   auto operator()(F && f, T const & t, Args && ... args) DECL_AND_RETURN( std::forward<F>(f)(args...));
 };

 template<typename F, typename T>
  auto apply_on_tuple (F && f, T const & t) DECL_AND_RETURN( apply_on_tuple_impl<std::tuple_size<T>::value-1,T>()(std::forward<F>(f),t));
}

#endif

