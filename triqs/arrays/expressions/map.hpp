
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

#ifndef TRIQS_ARRAYS_EXPRESSION_MAP_H
#define TRIQS_ARRAYS_EXPRESSION_MAP_H
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include "../impl/common.hpp"
#include "./function_objects.hpp"
/* Given a function f : arg_type -> result_type, map(f) is the function promoted to arrays
   map(f) : array<arg_type, N, Opt> --> array<result_type, N, Opt> 
   */
namespace triqs { namespace arrays { 
 
 template<class F> class map_impl  { 
  F f;
  
  template<typename A> 
   class map1_result : TRIQS_MODEL_CONCEPT(ImmutableArray) { 
    public:
     typedef typename F::result_type value_type;
     typedef typename A::domain_type domain_type;
     A const & a; F f;
     map1_result(F const & f_, A const & a_):a(a_),f(f_) {}
     domain_type domain() const { return a.domain(); } 
     value_type operator[] ( typename domain_type::index_value_type const & key) const { return f(a[key]); }
   };

  template<class A, class B> 
   class map2_result : TRIQS_MODEL_CONCEPT(ImmutableArray) { 
    public:
     typedef typename F::result_type value_type;
     typedef typename A::domain_type domain_type;
     A const & a; B const & b; F f;
     map2_result(F const & f_, A const & a_, B const & b_):a(a_),b(b_),f(f_) {
      if (a.domain() != b.domain()) TRIQS_RUNTIME_ERROR<<"map2 : domain mismatch";
     }
     domain_type domain() const { return a.domain(); } 
     value_type operator[] ( typename domain_type::index_value_type const & key) const { return f(a[key],b[key]); }
   };

  template <class T1, class V> struct argcheck : 
   boost::is_same< typename boost::remove_const<typename boost::remove_reference<T1>::type>::type , V>{};

  public : 
  map_impl(F const & f_):f(f_) {}

  template< class A >
   typename boost::enable_if_c< (F::arity==1), map1_result<A> >::type 
   operator()(A const & a) { 
    static_assert( (ImmutableArray<A>::value), "map1 : A does not have ImmutableArray");
    static_assert( (argcheck< typename F::arg1_type , typename A::value_type>::value), "type mismatch");
    return map1_result<A>(f,a);
   } 

  template< class A, class B >
   typename boost::enable_if_c< (F::arity==2), map2_result<A,B> >::type 
   operator()(A const & a, B const & b) { 
    static_assert( (ImmutableArray<A>::value), "map1 : A does not have ImmutableArray");
    static_assert( (ImmutableArray<B>::value), "map1 : B does not have ImmutableArray");
    static_assert( (argcheck< typename F::arg1_type , typename A::value_type>::value), "type mismatch");
    static_assert( (argcheck< typename F::arg2_type , typename A::value_type>::value), "type mismatch");
    static_assert( (boost::is_same<typename  A::domain_type, typename  B::domain_type>::value), "type mismatch");
    return map2_result<A,B>(f,a,b);
   } 
 };

 template<class F> map_impl<F> map (F const & f) { return map_impl<F>(f); }

 template<typename R, typename A1> 
  map_impl<function_object::from_regular_function< R (*)(A1) > >
  map (R (*f)(A1) ) { return map (function_object::make_from_regular_function(f));}

 template<typename R, typename A1, typename A2> 
  map_impl<function_object::from_regular_function< R (*)(A1,A2) > >
  map (R (*f)(A1,A2) ) { return map (function_object::make_from_regular_function(f));}

 template<class F> std::ostream & operator<<(std::ostream & out, map_impl<F> const & x){ return out<<"map("<<"F"<<")";}

 namespace result_of { 
  template<class F> struct map { typedef map_impl<F> type;};
  template<typename R, typename A1> 
   struct  map <R (*)(A1) > { typedef map_impl<function_object::from_regular_function< R (*)(A1) > >  type;};
  template<typename R, typename A1, typename A2> 
   struct  map <R (*)(A1,A2) > { typedef map_impl<function_object::from_regular_function< R (*)(A1,A2) > >  type;};
 }

}}//namespace triqs::arrays

#endif

