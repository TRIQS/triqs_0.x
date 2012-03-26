/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by O. Parcollet
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
#ifndef TRIQS_UTILITY_ALGEBRA_H
#define TRIQS_TAYLOR_EXPANSION_H 

#include <boost/utility/enable_if.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/proto/core.hpp>
#include <boost/proto/context.hpp>
#include <boost/proto/transform/arg.hpp>
#include <boost/proto/transform.hpp>
//#include <boost/preprocessor/repetition/repeat.hpp>
//#include <boost/preprocessor/repetition/repeat_from_to.hpp>
//#include <boost/preprocessor/facilities/intercept.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
//#include <boost/preprocessor/punctuation/comma.hpp>
//#include <boost/preprocessor/arithmetic/sub.hpp>
#include <complex>
//#include "triqs/utility/typeid_name.hpp"
#include <assert.h>

namespace triqs { namespace utility { namespace expressions { 

 namespace proto = boost::proto; namespace p_tag= proto::tag;

 template<typename T> struct is_in_ZRC : boost::is_arithmetic<T>  {};
 template<> struct is_in_ZRC<bool> : mpl::true_ {};
 template<typename T> struct is_in_ZRC<std::complex<T> > :  mpl::true_ {};

 // technical : evaluation of arithmetic operators
 template<typename TAG, typename A, typename B> struct _ops_;

 // the list of allowed operators in the grammar
#define BINARY_OP ((plus,+))((minus,-))((multiplies,*))((divides,/))

#define OP_NAME(elem) BOOST_PP_TUPLE_ELEM(2,0,elem)
#define OP_OP(elem)   BOOST_PP_TUPLE_ELEM(2,1,elem)

 // a trick to avoid putting T() in the typeof type deduction ! This code is NEVER used 
 template<class T> typename boost::unwrap_reference<T>::type pseudo_default_construct() { 
  typename boost::unwrap_reference<T>::type * x= NULL; assert(0); return *x; 
 }

#define AUX(r, data, elem) \
 template<typename A, typename B> struct _ops_<p_tag::OP_NAME(elem),A,B> {\
  typedef BOOST_TYPEOF_TPL( pseudo_default_construct<A>() OP_OP(elem) pseudo_default_construct<B>()) result_type;\
  static result_type invoke( A const &a , B const & b) { return a OP_OP(elem) b;} };
 BOOST_PP_SEQ_FOR_EACH(AUX, nil , BINARY_OP); 
#undef AUX

 template <typename T> std::ostream & formal_print(std::ostream & out, T const & x) { return out<<x;}

 template< 
  template<typename T> class leaf_identification >
  ,template<typename T> class scalar_identification >
  ,template<typename T> class wrap_scalar >
  ,template<typename ProtoTag, typename L, typename R> class wrap_binary_node >
  ,template<typename T> class wrap_negate >
  >
  struct grammar_generator {

   struct LeafGrammar :  proto::and_< proto::terminal<proto::_>, proto::if_<leaf_identification<proto::_value>()> > {}; 

   struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<scalar_identification<proto::_value>()> > {}; 

   //template <typename Expr> struct TheExpr;

   struct Grammar : 
    proto::or_<
    proto::when< ScalarGrammar,wrap_scalar<proto::_value>(proto::_value) >
    ,proto::when< LeafGrammar,proto::_value >
    ,proto::when< proto::plus <Grammar,Grammar>,      wrap_binary_node<p_tag::plus,_left,_right > (_left,_right) >
    ,proto::when< proto::minus <Grammar,Grammar>,     wrap_binary_node<p_tag::minus,_left,_right > (_left,_right)  >
    ,proto::when< proto::multiplies<Grammar,Grammar>, wrap_binary_node<p_tag::multiplies,_left,_right > (_left,_right)>
    ,proto::when< proto::divides<Grammar,Grammar>,    wrap_binary_node<p_tag::divides,_left,_right > (_left,_right)>
    ,proto::when< proto::negate<Grammar >,            wrap_negate <_left >(_left) >
    > {};

  };

 /* -------------------------------------------
  *   Print context
  * ------------------------------------------ */

 struct AlgebraPrintCtx : proto::callable_context< AlgebraPrintCtx const > {
  typedef std::ostream &result_type;
  result_type out;
  AlgebraPrintCtx(std::ostream & out_):out(out_) {}
  template <typename T>
   result_type operator ()(proto::tag::terminal, const T & A) const { return formal_print(out,A); }
  template<typename L, typename R>
   result_type operator ()(proto::tag::plus, L const &l, R const &r) const { return out << '(' << l << " + " << r << ')'; }
  template<typename L, typename R>
   result_type operator ()(proto::tag::minus, L const &l, R const &r) const { return out << '(' << l << " - " << r << ')'; }
  template<typename L, typename R>
   result_type operator ()(proto::tag::multiplies, L const &l, R const &r) const { return out << l << " * " << r; }
  template<typename L, typename R>
   result_type operator ()(proto::tag::divides, L const &l, R const &r) const { return out << l << " / " << r; }
 };


 // given an algebra A, build the algebra of A-valued functions.

 namespace algebra_function { 

  template<typename ProtoTag, typename L, typename R> struct wrap_binary_node  { 
   L const & l; R const & r;
   wrap_binary_node (L const & l_, R const & r_):l(l_),r(r_) {}
   template <typename T> struct call_rtype {
    typedef BOOST_TYPEOF_TPL (pseudo_default_construct<L>() (pseudo_default_construct<T>())) T1;
    typedef BOOST_TYPEOF_TPL (pseudo_default_construct<R>() (pseudo_default_construct<T>())) T2;
    typedef _ops_<ProtoTag, typename L:: template call_rtype<T>::type, typename R:: template call_rtype<R>::type > ops_type;
    typedef typename ops_type::result_type type;
   };
   template<typename T> typename call_rtype<T>::type  operator() (T const & arg) const { return  call_rtype<T>::ops_type::invoke(l(arg),r(arg));}
  }; 

  template<typename S> struct wrap_scalar {
   typedef T result_type;
   S s; 
   wrap_scalar( S const & x) : s(x) {}
   template<typename T> result_type operator() (T) const { return s;}
  };

  template<typename L> struct wrap_negate  { 
   L const & l; 
   wrap_negate (L const & l_):l(l_) {} 
   template <typename T> struct call_rtype {
    typedef BOOST_TYPEOF_TPL ( (- pseudo_default_construct<L>() (pseudo_default_construct<T>()))) type;
   };
   template<typename T> typename call_rtype<T>::type operator() (T const & arg) const { return  (- l(arg));} 
  }; 

  template< 
  template<typename T> class leaf_identification >
  ,template<typename T> class algebra_element_identification >
  >
  struct function_grammar_generator: grammar_generator< leaf_identification, algebra_element_identification, wrap_scalar, wrap_binary_node,wrap_negate> {};

 }

  //static_assert( !(boost::is_same<ProtoTag,p_tag::greater >::value), "oops");
 //greater_equal

}}}

#undef OP_NAME
#undef OP_OP
#undef BINARY_OP

#endif

