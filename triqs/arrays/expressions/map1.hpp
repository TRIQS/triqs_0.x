
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

#ifndef EXPRESSION_MAP_H
#define EXPRESSION_MAP_H 

#include <boost/utility/enable_if.hpp>
#include "../matrix.hpp"
#include "../array.hpp"
#include "./make_function_macros.hpp"

namespace triqs { namespace arrays { 
 namespace Expression { namespace details { template<class F, class A> struct map_impl1; template<class F, class A, class B> struct map_impl2; }} // impl below

 namespace result_of { 
  template<class F, typename A, typename B=void> struct map {  
   //static_assert( (is_matrix_or_view<A>::value), "map : the argument must be a matrix"); 
   typedef Expression::details::map_impl2<F,A,B> type;}; 

  template<class F, typename A> struct map<F,A,void> {  
   //static_assert( (is_matrix_or_view<A>::value), "map : the argument must be a matrix"); 
   typedef Expression::details::map_impl1<F,A> type;}; 
 }

 namespace Expression {

  template<class F, typename A> 
   typename result_of::map<F,A,void>::type map (F const & f, A const & a) { return details::map_impl1<F,A>(f,a); }

  template<class F, typename A,typename B> 
   typename result_of::map<F,A,B>::type map (F const & f, A const & a, B const & b) { return details::map_impl2<F,A,B>(f,a,b); }

#define ARRAY_EXPRESSION_MAP_CALLABLE_OBJECT1(Name,OBJ) \
  template<class T> typename boost::enable_if<has_immutable_array_interface<T>,typename result_of::map<OBJ,T>::type >::type\
  Name(T const & x) { return Expression::map( OBJ(), x); }

#define ARRAY_EXPRESSION_MAP_FUNCTION1(Name,F) ARRAY_EXPRESSION_MAKE_FUNCTION1(Name,F); ARRAY_EXPRESSION_MAP_CALLABLE_OBJECT1(Name,Name##_impl) 

#define ARRAY_EXPRESSION_MAP_CALLABLE_OBJECT2(Name,OBJ) \
  template<class T1, class T2> typename boost::enable_if<has_immutable_array_interface2<T1,T2> , typename result_of::map<OBJ,T1,T2>::type >::type\
  Name(T1 const & x, T2 const & y) { return Expression::map( OBJ(), x,y); }

#define ARRAY_EXPRESSION_MAP_FUNCTION2(Name,F) ARRAY_EXPRESSION_MAKE_FUNCTION2(Name,F); ARRAY_EXPRESSION_MAP_CALLABLE_OBJECT2(Name,Name##_impl) 

  namespace details {//------------- IMPLEMENTATION -----------------------------------

   template<class F, typename A> 
    class map_impl1 : Tag::expression_terminal, Tag::has_immutable_array_interface, Tag::expression {
     typedef typename A::domain_type::index_value_type index_value_type;

     public:
     typedef typename A::value_type value_type;
     typedef typename A::domain_type domain_type;
     A const & a; F f;

     map_impl1(F const & f_, A const & a_):a(a_),f(f_) {}
     domain_type domain() const { return a.domain(); } 
     value_type operator[] (index_value_type const & key) const { return f(a[key]); }
    };

   template<class F, typename A, typename B> 
    class map_impl2 : Tag::expression_terminal, Tag::has_immutable_array_interface, Tag::expression {
     typedef typename A::domain_type::index_value_type index_value_type;

     public:
     typedef typename A::value_type value_type;
     typedef typename A::domain_type domain_type;
     A const & a; B const & b; F f;

     map_impl2(F const & f_, A const & a_, B const & b_):a(a_),b(b_),f(f_) {}
     domain_type domain() const { return a.domain(); } 
     value_type operator[] (index_value_type const & key) const { return f(a[key],b[key]); }
    };

   template<class F,typename A> 
    std::ostream & operator<<(std::ostream & out, map_impl1<F,A> const & x){ return out<<"map(F"<<","<<x.a<<")";}
   template<class F,typename A,typename B> 
    std::ostream & operator<<(std::ostream & out, map_impl2<F,A,B> const & x){ return out<<"map(F"<<","<<x.a<<","<<x.b<<")";}

  } // details
 } // Expression
} // namespace triqs_arrays 
#endif

