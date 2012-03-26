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
#define TRIQS_UTILITY_ALGEBRA_H 

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

namespace triqs { namespace utility { namespace proto { 

  namespace mpl = boost::mpl; namespace proto = boost::proto; namespace p_tag= proto::tag;

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

 namespace algebra { 

 template< typename A, template<typename Expr> class The_Expr> 
  struct grammar_generator {

   struct LeafGrammar   : proto::and_< proto::terminal<proto::_>, proto::if_<typename A::template is_element<proto::_value>()> > {}; 
   struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<typename A::template is_scalar<proto::_value>()> > {}; 
   struct Grammar : 
    proto::or_<
    proto::when< ScalarGrammar,                       typename A::template wrap_scalar<proto::_value>(proto::_value) >
    ,proto::when< LeafGrammar,                        proto::_value >
    ,proto::when< proto::plus <Grammar,Grammar>,      typename A::template wrap_binary_node<p_tag::plus,proto::_left,proto::_right > (proto::_left,proto::_right) >
    ,proto::when< proto::minus <Grammar,Grammar>,     typename A::template wrap_binary_node<p_tag::minus,proto::_left,proto::_right > (proto::_left,proto::_right)  >
    ,proto::when< proto::multiplies<Grammar,Grammar>, typename A::template wrap_binary_node<p_tag::multiplies,proto::_left,proto::_right > (proto::_left,proto::_right)>
    ,proto::when< proto::divides<Grammar,Grammar>,    typename A::template wrap_binary_node<p_tag::divides,proto::_left,proto::_right > (proto::_left,proto::_right)>
    ,proto::when< proto::negate<Grammar >,            typename A::template wrap_negate <proto::_left >(proto::_left) >
    > {};

   typedef Grammar type;
  };

 template< typename Grammar, template<typename Expr> class The_Expr, bool WithCopy> struct domain;

 template< typename Grammar, template<typename Expr> class The_Expr> struct domain<Grammar,The_Expr,false> : 
  proto::domain<proto::generator<The_Expr>, Grammar> { };
 template< typename Grammar, template<typename Expr> class The_Expr> struct domain<Grammar,The_Expr,true> : 
  proto::domain<proto::generator<The_Expr>, Grammar> {
  //template< typename T > struct as_child : proto_base_domain::as_expr< T > {};
 };

 }

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


 }}}

#undef OP_NAME
#undef OP_OP
#undef BINARY_OP

#endif

