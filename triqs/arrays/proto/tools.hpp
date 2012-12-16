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
#ifndef TRIQS_ARRAYS_EXPRESSION_TOOLS_H
#define TRIQS_ARRAYS_EXPRESSION_TOOLS_H
#include "../array.hpp"
namespace triqs { namespace arrays {

 namespace tags { struct plus{}; struct minus{}; struct multiplies{}; struct divides{}; struct is_scalar{}; }

// The basic operations put in a template.... 
 template<typename Tag> struct operation;
 template<> struct operation<tags::plus> { 
  template<typename L, typename R> auto operator()(L const & l, R const & r) const -> decltype(l+r) { return l+r;} 
  static const char name = '+';
 };
 template<> struct operation<tags::minus> { 
  template<typename L, typename R> auto operator()(L const & l, R const & r) const -> decltype(l-r) { return l-r;} 
  static const char name = '-';
 };
 template<> struct operation<tags::multiplies> { 
  template<typename L, typename R> auto operator()(L const & l, R const & r) const -> decltype(l*r) { return l*r;} 
  static const char name = '*';
 };
 template<> struct operation<tags::divides> { 
  template<typename L, typename R> auto operator()(L const & l, R const & r) const -> decltype(l/r) { return l/r;} 
  static const char name = '/';
 };

// The scalar ...
 template<typename T> struct is_in_ZRC : boost::is_arithmetic<T>  {};
 template<> struct is_in_ZRC<bool> : mpl::true_ {};
 template<typename T> struct is_in_ZRC<std::complex<T> > :  mpl::true_ {};

 // Wrapping the scalar in a little struct to recognize it
 template<typename S> struct array_scalar_wrap : tags::is_scalar {
  typedef S value_type; 
  S s; array_scalar_wrap(S const &s_):s(s_){} 
  template<typename KeyType> value_type operator[](KeyType&& key) const {return s;}
  friend std::ostream &operator <<(std::ostream &sout, array_scalar_wrap const &expr){return sout << expr.s; }
 };

 // get the rank of something ....
 template<typename T> struct get_rank { static constexpr int value = T::domain_type::rank;};
 template<typename S> struct get_rank<array_scalar_wrap<S>> { static constexpr int value =0;};
 
 //
 template<typename T> struct keeper_type : boost::mpl::if_<is_in_ZRC<T>, array_scalar_wrap<T>, typename view_type_if_exists_else_type<T>::type> {};

 // Combine the two domains of LHS and RHS : need to specialize where there is a scalar
 struct combine_domain {
  template<typename L, typename R> 
   auto operator() (L const & l, R const & r) const -> decltype(l.domain()) { 
    if (l.domain().lengths() != r.domain().lengths()) TRIQS_RUNTIME_ERROR << "Domain size mismatch : "<< l.domain().lengths()<<" vs" <<r.domain().lengths();
    return l.domain();
   }
  template<typename S, typename R> auto operator() (array_scalar_wrap<S> const & w, R const & r) const -> decltype(r.domain()) { return r.domain();}
  template<typename S, typename L> auto operator() (L const & l, array_scalar_wrap<S> const & w) const -> decltype(l.domain()) { return l.domain();}
 };
 
}}//namespace triqs::arrays
#endif
