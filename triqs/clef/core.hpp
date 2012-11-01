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
#ifndef TRIQS_CLEF_CORE_H
#define TRIQS_CLEF_CORE_H

#define TRIQS_HAS_CLEF_EXPRESSIONS
#define TRIQS_CLEF_MAXNARGS 8
#define TRIQS_CLEF_MAXNARGS_CALLABLE TRIQS_CLEF_MAXNARGS

#include <triqs/utility/first_include.hpp>

#ifdef TRIQS_ARRAYS_ALREADY_INCLUDED
#error "If you use triqs::clef and triqs::arrays, you MUST include clef first"
#endif

#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/filter_view.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/proto/core.hpp>
#include <boost/proto/context.hpp>
#include <boost/proto/transform/arg.hpp>
#include <boost/proto/transform.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/map_fwd.hpp>
#include <boost/fusion/include/has_key.hpp>
#include <boost/fusion/include/make_map.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_trailing_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_shifted_params.hpp>
#include <boost/preprocessor/repetition/enum_params_with_a_default.hpp>
#include <boost/preprocessor/facilities/intercept.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <complex>
//#include "triqs/utility/typeid_name.hpp"
#include <assert.h>
//#include <iostream> // debug only 

#include <algorithm> // debug print of ph only
#include <iterator> // debug print of ph only
#include <vector> 

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#define USE_STATIC_ASSERT
#define __USE_C11__
#endif

#ifndef USE_STATIC_ASSERT
#include "boost/static_assert.hpp"
#define static_assert(X,MESS) BOOST_STATIC_ASSERT((X)) 
#endif

// the list of allowed operators in the grammar
#define ALLOWED_BINARY_OP \
 ((plus,+))((minus,-))((multiplies,*))((divides,/))((greater,>))((less,<))((greater_equal,>=))((less_equal,<=))((equal_to,==))

#define OP_NAME(elem) BOOST_PP_TUPLE_ELEM(2,0,elem)
#define OP_OP(elem)   BOOST_PP_TUPLE_ELEM(2,1,elem)

#define typename_A(N)              BOOST_PP_ENUM_PARAMS(N,typename A) 
#define _A_const_ref(N)            BOOST_PP_ENUM_BINARY_PARAMS(N,A,const & BOOST_PP_INTERCEPT)
#define _A_const_ref_a(N)          BOOST_PP_ENUM_BINARY_PARAMS(N,A,const & a)
#define _A_const_ref_a_tr(N)       BOOST_PP_ENUM_TRAILING_BINARY_PARAMS(N,A,const & a)
#define _A_ref(N)                  BOOST_PP_ENUM_BINARY_PARAMS(N,A, & BOOST_PP_INTERCEPT)
#define _A_(N)                     BOOST_PP_ENUM_PARAMS(N,A)
#define _A_void(N)                 BOOST_PP_ENUM_PARAMS(N,A) BOOST_PP_ENUM_TRAILING_PARAMS(BOOST_PP_SUB(TRIQS_CLEF_MAXNARGS, N), void BOOST_PP_INTERCEPT)
#define _a_(N)                     BOOST_PP_ENUM_PARAMS(N,a)

namespace triqs { namespace clef { 

 namespace mpl = boost::mpl; namespace bf = boost::fusion; namespace proto = boost::proto; namespace p_tag= proto::tag;
 using boost::enable_if; using boost::disable_if;
 using proto::if_else; // adding this function to this namespace
 /* ---------------------------------------------------------------------------------------------------
  * Define the main expression template FntExpr, the domain and Grammar (implemented below)
  * NB : it modifies the PROTO Domain to make *COPIES* of all objects.
  * cf http://www.boost.org/doc/libs/1_49_0/doc/html/proto/users_guide.html#boost_proto.users_guide.front_end.customizing_expressions_in_your_domain.per_domain_as_child
  * Objects which have a view type are however NOT copied, the copy is replaced by a VIEW.
  *   --------------------------------------------------------------------------------------------------- */
 template <typename Expr> struct FntExpr;

 template<typename T, typename Void =void> struct view_type_if_exists_else_type {typedef T type;}; 
 template<typename T> struct view_type_if_exists_else_type<T, typename T::has_view_type_tag> {typedef const typename T::view_type type;}; 

 struct FntGrammar;
 struct FntDomain : proto::domain<proto::generator<FntExpr>, FntGrammar> {
  //template< typename T > struct as_child : proto_base_domain::as_expr< T > {};
  template< typename T > struct as_child : proto_base_domain::as_expr<typename boost::add_const<typename view_type_if_exists_else_type <T>::type >::type > {};
 };

 /* ---------------------------------------------------------------------------------------------------
  *  Placeholder and corresponding traits
  *  --------------------------------------------------------------------------------------------------- */
 struct generic_domain_type {};

 template<int N> struct placeholder { 
  template <typename RHS> bf::pair<placeholder,RHS> operator = (RHS rhs) { return bf::make_pair<placeholder>(boost::cref(rhs)); } 
 };

 template<typename T> struct is_placeholder : mpl::false_  {}; 
 template<int N> struct is_placeholder <placeholder <N> > : mpl::true_  {}; 

 template<typename PH, typename BfMapType> struct ph_in_map :  mpl::has_key<BfMapType, PH>::type {};

 template<typename T> struct is_proto_terminal : is_placeholder<T> {}; // only placeholders are automatically terminals
 template<typename T> struct is_proto_terminal<boost::reference_wrapper<T> > : is_proto_terminal<T> {};

 /* ---------------------------------------------------------------------------------------------------
  *  Detect is an object is a lazy expression
  *  --------------------------------------------------------------------------------------------------- */

#define TRIQS_CLEF_IS_EXPRESSION() typedef void __triqs_nvl_is_expression_tag;

 template <typename T, typename Void =void> struct is_custom_lazy_expression : mpl::false_{};
 template <typename T> struct is_custom_lazy_expression<T, typename T::__triqs_nvl_is_expression_tag> : mpl::true_{};

 template <typename T> struct is_lazy : mpl::or_<is_custom_lazy_expression<T>, is_placeholder< T> > {};
 template <> struct is_lazy<void> : mpl::false_ {}; 
 template <typename T> struct is_lazy< FntExpr<T> > : mpl::true_ {};
 template <typename T> struct is_lazy< T&& > : is_lazy<T> {};
 template <typename T> struct is_lazy< T& > : is_lazy<T> {};

 template<typename T> struct is_not_lazy : mpl::not_< is_lazy <T> > {}; 

 //#if __has_feature(cxx_variadic_templates) // clang only 
#ifdef __USE_C11__
 template<typename... Args> struct one_is_lazy;
 template<typename T, typename... Args> struct one_is_lazy<T, Args...> : 
  boost::mpl::or_ < triqs::clef::is_lazy<T>, one_is_lazy<Args...> > {};
 template<> struct one_is_lazy<> : boost::mpl::false_ {};
#else
 // check that TRIQS_CLEF_MAXNARGS is < 10
 template< BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(10, typename Arg, void ) >  struct one_is_lazy : 
  boost::mpl::or_< 
  boost::mpl::or_< is_lazy<Arg0>, is_lazy<Arg1>, is_lazy<Arg2>, is_lazy<Arg3>, is_lazy<Arg4> >,
  boost::mpl::or_< is_lazy<Arg5>, is_lazy<Arg6>, is_lazy<Arg7>, is_lazy<Arg8>, is_lazy<Arg9> > > {};
#endif

 /* *****************************************************************************************
  * Grammar 
  * *****************************************************************************************/

 struct NonLazyLeafGrammar       : proto::and_< proto::terminal<proto::_>, proto::if_<is_not_lazy<proto::_value>()> > {}; 
 struct PlaceholderLeafGrammar   : proto::and_< proto::terminal<proto::_>, proto::if_<is_placeholder<proto::_value>()> > {}; 
 struct CustomExprLeafGrammar      : proto::and_< proto::terminal<proto::_>, proto::if_<is_custom_lazy_expression<proto::_value>()> > {}; 

 struct CallableGrammar          : proto::or_ < NonLazyLeafGrammar , FntGrammar> {}; 
 struct SubscriptableGrammar     : proto::or_ < NonLazyLeafGrammar , FntGrammar> {}; 

 struct FntGrammarCases { template <typename TAG> struct case_: proto::not_<proto::_>{}; };
 template<> struct FntGrammarCases::case_<proto::tag::terminal>:proto::or_< NonLazyLeafGrammar, PlaceholderLeafGrammar, CustomExprLeafGrammar > {};

 // Matching all the binary operations allowed in the grammar
 // We do NOT match NonLazyLeafGrammar OP NonLazyLeafGrammar, it has to be done by the operator of the real object ....
#define ONETYPE(r, data, elem)  \
 template<> struct FntGrammarCases::case_<proto::tag::OP_NAME(elem)> : \
 proto::and_< proto::OP_NAME(elem) <FntGrammar,FntGrammar>, proto::not_<proto::OP_NAME(elem) <NonLazyLeafGrammar,NonLazyLeafGrammar> > > {};
 BOOST_PP_SEQ_FOR_EACH(ONETYPE, nil , ALLOWED_BINARY_OP); 
#undef ONETYPE

 template<> struct FntGrammarCases::case_<proto::tag::function> :proto::function<CallableGrammar, proto::vararg<FntGrammar> >{};
 template<> struct FntGrammarCases::case_<proto::tag::subscript>: proto::subscript <SubscriptableGrammar,FntGrammar> {}; 
 template<> struct FntGrammarCases::case_<proto::tag::if_else_>: proto::if_else_ <FntGrammar,FntGrammar,FntGrammar> {};  
 struct FntGrammar : proto::switch_<FntGrammarCases> {};
 
 /* ---------------------------------------------------------------------------------------------------
  * a metafunction to extract the type of a bf::map with a default if it is not present.
  * --------------------------------------------------------------------------------------------------- */
 template<typename BfMapType, typename Key, typename Default, typename Enable=void>
  struct bf_result_of_at_key_w_default { typedef Default type;};

 template<typename BfMapType, typename Key, typename Default>
  struct bf_result_of_at_key_w_default<BfMapType,Key,Default, typename enable_if< mpl::has_key<BfMapType,Key> >::type > : 
  bf::result_of::at_key<BfMapType const ,Key> {};

 /* ---------------------------------------------------------------------------------------------------
  * eval_t : first the node transformations 
  * --------------------------------------------------------------------------------------------------- */

 struct eval_node_t { // any leaf which is not a ph
  BOOST_PROTO_CALLABLE();
  template<typename Sig> struct result;
  template<typename This, typename T> struct result<This(T)> { typedef typename boost::remove_const<typename boost::remove_reference<T>::type>::type type;};
  template<typename This, typename T> struct result<This(boost::reference_wrapper<T>)> { typedef T const & type;};
  template<typename This, typename T> struct result<This(boost::reference_wrapper<T> const &)> { typedef T const & type;};

  template<typename T> typename result<eval_node_t(T)>::type operator ()(T const & x) const { return x;}

  template<typename T> typename result<eval_node_t(boost::reference_wrapper<T> const &)>::type 
   operator ()(boost::reference_wrapper<T> const & x) const { return boost::unwrap_ref(x); } 
 };

 struct eval_node_ph_t { // ph leaf
  BOOST_PROTO_CALLABLE();
  template<typename Sig> struct result;
  template<typename This, typename PH, typename BfMapType> struct result<This(PH,BfMapType)> {
   typedef typename boost::remove_const<typename boost::remove_reference<BfMapType>::type>::type B;
   typedef typename boost::remove_const<typename boost::remove_reference<PH>::type>::type P;
   typedef typename bf_result_of_at_key_w_default< B, P , typename  proto::result_of::make_expr<proto::tag::terminal, FntDomain, P const >::type >::type type; 
  };

  template<typename PH, typename BfMapType>
   typename enable_if< ph_in_map<PH,BfMapType>, typename result<eval_node_ph_t(PH,BfMapType)>::type >::type  
   operator ()(PH const & , BfMapType const & ph_value_dict) const { return bf::at_key<PH>(ph_value_dict); }

  template<typename PH, typename BfMapType>
   typename disable_if< ph_in_map<PH,BfMapType>,  typename result<eval_node_ph_t(PH,BfMapType)>::type >::type  
   operator ()(PH const &, BfMapType const & ph_value_dict) const { return proto::make_expr<proto::tag::terminal,FntDomain>(PH()); }
 };

 struct eval_custom_expr_t { // custom expression
  BOOST_PROTO_CALLABLE();
  template<typename Sig> struct result;
  template<typename This, typename Expr, typename BfMapType> struct result<This(Expr,BfMapType)> {
   typedef typename boost::remove_const<typename boost::remove_reference<BfMapType>::type>::type B;
   typedef typename boost::remove_const<typename boost::remove_reference<Expr>::type>::type E;
   typedef typename E::template call_rtype<B>::type type; 
  };
  template<typename Expr, typename BfMapType> typename result<eval_node_ph_t(Expr,BfMapType)>::type
   operator ()(Expr const & expr, BfMapType const & ph_value_dict) const { return expr.eval_on_map(ph_value_dict); }
 };

 // a trick to avoid putting T() in the typeof type deduction ! This code is NEVER used. It is a hack to avoid decltype....
 template<class T> typename boost::remove_reference<typename boost::unwrap_reference<T>::type>::type * pseudo_default_construct() { 
  typename boost::remove_reference<typename boost::unwrap_reference<T>::type>::type * x= NULL; assert(0); return x; 
 }

 // transform to eval the function call 
#define AUX6(z,p,unused) (*(pseudo_default_construct<typename boost::remove_reference<A##p>::type>()))
#define AUX(z, NN, unused) \
 struct call_fnt_t##NN {\
  BOOST_PROTO_CALLABLE();\
  template<typename Sig> struct result;\
  template<typename This, typename F BOOST_PP_ENUM_TRAILING_PARAMS(NN, typename A)> struct result<This(F BOOST_PP_ENUM_TRAILING_PARAMS(NN,A))> {\
   typedef BOOST_TYPEOF_TPL ( (*(pseudo_default_construct<const F>())) (BOOST_PP_ENUM(NN,AUX6,nil))) type;\
  };\
  template<typename F BOOST_PP_ENUM_TRAILING_PARAMS(NN, typename A)>\
  typename result<call_fnt_t##NN(F BOOST_PP_ENUM_TRAILING_PARAMS(NN,A))>::type\
  operator()(const F & f BOOST_PP_ENUM_TRAILING_BINARY_PARAMS(NN,A, const &a)) const { return f(BOOST_PP_ENUM_PARAMS(NN,a));}\
 };
 BOOST_PP_REPEAT(BOOST_PP_INC(TRIQS_CLEF_MAXNARGS_CALLABLE), AUX, nil)
#undef AUX
#undef AUX6

  /* ---------------------------------------------------------------------------------------------------
   * eval_t : the proto transform that makes the evaluations. Almost trivial, except for fnt call
   * boost proto has trouble finding hte result of a call (without implementing ResultOf concept in the user object
   * which I do not want to impose...)
   * --------------------------------------------------------------------------------------------------- */

  struct eval_t; 
 struct eval_t_cases { template <typename TAG> struct case_: proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::_default<eval_t > > {}; };

 template<> struct eval_t_cases::case_<proto::tag::terminal>:
  proto::or_<
  proto::when<NonLazyLeafGrammar, eval_node_t(proto::_value)>
  ,proto::when<PlaceholderLeafGrammar, eval_node_ph_t(proto::_value, proto::_state) >
  ,proto::when<CustomExprLeafGrammar, eval_custom_expr_t(proto::_value, proto::_state) >
  >{};

 template<> struct eval_t_cases::case_<proto::tag::function>:
  proto::or_<
#define AUX2(z,p,unused) eval_t(proto::_child##p)
#define AUX(z, NN, unused)  \
  proto::when< proto::function <CallableGrammar BOOST_PP_ENUM_TRAILING_PARAMS(NN,FntGrammar BOOST_PP_INTERCEPT)>, \
  call_fnt_t##NN (BOOST_PP_ENUM(BOOST_PP_INC(NN), AUX2, nil )) >
  BOOST_PP_ENUM(BOOST_PP_INC(TRIQS_CLEF_MAXNARGS_CALLABLE) , AUX, nil)
#undef AUX
#undef AUX2
  >{};

 struct eval_t : proto::switch_<eval_t_cases> {};

 /* ---------------------------------------------------------------------------------------------------
  * eval
  * --------------------------------------------------------------------------------------------------- */
 namespace result_of {  template < typename C, typename BfMapType> struct eval_on_map : C::template call_rtype<BfMapType > {}; }

 template<typename C, typename BfMapType> 
  typename result_of::eval_on_map<C,BfMapType >::type 
  eval_on_map(C const & expr, BfMapType const & bfm ) { 
   static_assert( is_lazy<C>::value, "Error : evaluation can only be called on expressions");
   return expr.eval_on_map(bfm) ; 
  }

#define AUX1(z,i,unused) bf::pair<X##i,V##i>
#define AUX2(z,i,unused) bf::pair<X##i,V##i> const & p##i
#define IMPL(z,NN,unused)  \
 template<typename C, BOOST_PP_ENUM_PARAMS(NN, typename X),BOOST_PP_ENUM_PARAMS(NN, typename V) > \
 typename result_of::eval_on_map<C,bf::map<BOOST_PP_ENUM(NN, AUX1,nil)> >::type \
 eval(C const & expr, BOOST_PP_ENUM(NN, AUX2,nil) ) {\
  return expr.eval_on_map(bf::map<BOOST_PP_ENUM(NN, AUX1,nil) > (BOOST_PP_ENUM_PARAMS(NN,p))) ; }
 BOOST_PP_REPEAT_FROM_TO (1, BOOST_PP_INC(TRIQS_CLEF_MAXNARGS) , IMPL , nil)
#undef AUX1
#undef AUX2
#undef IMPL

  /* *****************************************************************************************
   * A transform to accumulate the placeholders in an expression
   * *****************************************************************************************/

  struct ph_accu_t { // to handle the PH leaf 
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename PH, typename BFL> struct result<This(PH,BFL)> { 
    typedef typename boost::remove_const<typename boost::remove_reference<PH>::type>::type P;
    typedef typename boost::remove_const<typename boost::remove_reference<BFL>::type>::type B;
    typedef typename mpl::insert<B,P >::type type;
   };
   template<typename PH, typename BFL> result<ph_accu_t(PH,BFL)>  
    operator ()(const PH & ph, BFL const & bfl) const { return result<ph_accu_t(PH,BFL)>(); }
  };

 struct ph_set_custom_expr_t { // to handle any "custom" expression
  BOOST_PROTO_CALLABLE();
  template<typename Sig> struct result;
  template<typename This, typename Expr, typename BFL> struct result<This(Expr,BFL)> {
   typedef typename boost::remove_const<typename boost::remove_reference<BFL>::type>::type B;
   typedef typename boost::remove_const<typename boost::remove_reference<Expr>::type>::type E;
   typedef typename mpl::fold< typename E::ph_set , B, mpl::insert<mpl::_1,mpl::_2> >::type type;
  };
  template<typename Expr, typename BFL> typename result<ph_set_custom_expr_t(Expr,BFL)>::type
   operator ()(Expr const & expr, BFL const & bfl) const { typename result<ph_set_custom_expr_t(Expr,BFL)>::type();}
 };

 struct ph_accumulator_t : 
  proto::or_<
  proto::when< NonLazyLeafGrammar, proto::_state>
  ,proto::when< PlaceholderLeafGrammar, ph_accu_t(proto::_value, proto::_state)  > 
  ,proto::when<CustomExprLeafGrammar, ph_set_custom_expr_t(proto::_value, proto::_state) >
  ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::fold<proto::_, proto::_state, ph_accumulator_t >() >
  >{};

 /* ---------------------------------------------------------------------------------------------------
  *  Handle the ph_set, the set of placeholder of an expression
  *  --------------------------------------------------------------------------------------------------- */

 template<typename X, bool emp = is_custom_lazy_expression<X>::value > struct get_ph_set{ typedef  mpl::set<> type;};
 template<typename X> struct get_ph_set<X,true> { typedef typename X::ph_set type;};
 template<typename X> struct get_ph_set<FntExpr<X>,false> { 
  static_assert( (proto::matches<FntExpr<X>,ph_accumulator_t>::value), "CLEF : internal error get_ph_set");
  typedef typename boost::result_of< ph_accumulator_t(X,mpl::set<>) >::type type;
 };
 template<typename X> struct get_ph_set<const FntExpr<X>,false> : get_ph_set<FntExpr<X>,false> {};

 template< typename T> struct  has_no_ph : mpl::empty< typename get_ph_set<T>::type > {};

 /* ---------------------------------------------------------------------------------------------------
  *   Print context : to print nicely the expression
  * --------------------------------------------------------------------------------------------------- */

 template <typename T> typename enable_if<boost::is_arithmetic<T>, std::ostream>::type & 
  triqs_nvl_formal_print(std::ostream & out, T const & x) { return out<<x;}

 template <typename T> std::ostream & triqs_nvl_formal_print(std::ostream & out, boost::reference_wrapper<T> const & x) 
 { return triqs_nvl_formal_print(out,boost::unwrap_ref(x));}

 template <int N> std::ostream& triqs_nvl_formal_print(std::ostream & out, placeholder<N> const &x){return out<< "x"<<N;}

 struct PrintCtx : proto::callable_context< PrintCtx const > {
  typedef std::ostream &result_type;
  result_type out;
  PrintCtx(std::ostream & out_):out(out_) {}
  template <typename T>
   inline result_type operator ()(proto::tag::terminal, const T & A) const { return triqs_nvl_formal_print(out,A); }

#define AUX(rrr, data, elem) template<typename L, typename R>\
  result_type operator ()(proto::tag::OP_NAME(elem) , L const &l, R const &r) const { \
   return out << '(' << l << " "<< BOOST_PP_STRINGIZE(OP_OP(elem)) <<" " << r << ')'; }
  BOOST_PP_SEQ_FOR_EACH(AUX, nil , ALLOWED_BINARY_OP); 
#undef AUX

#define AUX1(z, p, unused) a##p << "," <<
#define AUX(z, NN, unused)  \
  template<typename F BOOST_PP_ENUM_TRAILING_PARAMS(NN,typename A) >\
  result_type operator ()(proto::tag::function, F const &f BOOST_PP_ENUM_TRAILING_BINARY_PARAMS(NN,A,const & a) ) const {\
   return out << f << '(' << BOOST_PP_REPEAT(NN,AUX1,nil) ')'; }
  BOOST_PP_REPEAT(BOOST_PP_INC(TRIQS_CLEF_MAXNARGS_CALLABLE) , AUX, nil)
#undef AUX
#undef AUX1

   template <typename F, typename T>
   result_type operator ()(proto::tag::subscript, F const & f, const T & A) const { return out<< f<< '[' <<A <<']'; }
 };

 /* ---------------------------------------------------------------------------------------------------
  * Debug code : print the list of placeholders, sorted alphabetically
  * --------------------------------------------------------------------------------------------------- */

 struct to_vec_ {
  std::vector<std::string> & v; to_vec_ (std::vector<std::string> & V_):v(V_){}
  template<typename T> void operator()(T) { std::stringstream out; triqs_nvl_formal_print( out, T()); v.push_back(out.str()); }
 };

 template <typename ph_set_type> std::string placeholder_list_as_string1() {
  std::vector<std::string> v; mpl::for_each<ph_set_type>(to_vec_(v)); 
  std::sort(v.begin(), v.end()); 
  std::stringstream fs; std::ostream_iterator<std::string> out_it (fs,",");
  std::copy(v.begin(), v.end(), out_it);
  return fs.str();
 }

 template <typename Expr> std::string placeholder_list_as_string( Expr expr) {
  return placeholder_list_as_string1<typename get_ph_set<Expr>::type>();
 }

 /* ---------------------------------------------------------------------------------------------------
  *  Implementation of the expression type FntExpr
  *  1) Expression can be evaluated using  : 
  *
  *    - eval_on_map (M) where M is a BfMapType i.e. a boost::fusion Associative Sequence
  *    - operator()  boost::fusion::pair objects
  *   2) FntExpr overload =  : expr = rhs, with two forms : 
  *      - A(x,y,...) = rhs where A is a assignable object and x,y are placeholder
  *      - x_ = value where x_ is a place holder : this return a boost::fusion::pair < PH_Type, value>
  *        this implies that expr( x_ = value will work as expected....)
  *        Note that value can be anything...
  * --------------------------------------------------------------------------------------------------- */
 // delegate the = operator later
 template <typename Expr, typename T, typename Enable=void > struct FntExpr_assign_delegation;

 template<typename Expr> struct FntExpr : proto::extends<Expr, FntExpr<Expr>, FntDomain>{
  FntExpr( Expr const & expr = Expr() ) : proto::extends<Expr, FntExpr<Expr>, FntDomain> ( expr ) {}

  template<typename BfMapType> struct call_rtype : boost::result_of<eval_t(FntExpr, BfMapType)> {};
  template<typename BfMapType> typename call_rtype<BfMapType>::type  eval_on_map(BfMapType const & bfm) const {return eval_t()(*this, bfm);}

  // I need to overload the [] and the () since the default of proto uses a ref for the _left... why is this ?
  // otherwise, it leads to a bug when chaining the calls like in f(z_)(x_,y_)
  template< typename Arg > 
   typename proto::result_of::make_expr< p_tag::subscript, triqs::clef::FntDomain ,FntExpr, Arg const &>::type  const 
   operator[]( Arg const & arg ) const { return proto::make_expr<p_tag::subscript, triqs::clef::FntDomain>( *this , arg);}

#define AUX(z,N, unused) \
  template< typename_A(N) > \
  typename proto::result_of::make_expr< p_tag::function, triqs::clef::FntDomain ,FntExpr,_A_const_ref(N)>::type const \
  operator()( _A_const_ref_a(N) ) const { return proto::make_expr<p_tag::function, triqs::clef::FntDomain>( *this, _a_(N) );}
  BOOST_PP_REPEAT_FROM_TO(1,BOOST_PP_INC(TRIQS_CLEF_MAXNARGS) , AUX, nil);
#undef AUX

  template<typename T> 
   typename FntExpr_assign_delegation<Expr ,T>::result_type
   operator=(T x) const { return FntExpr_assign_delegation<Expr ,T>::invoke(*this,x);}

  friend std::ostream &operator <<(std::ostream &sout, FntExpr<Expr> const &expr) { return proto::eval(expr, PrintCtx (sout)); }
  friend std::ostream & triqs_nvl_formal_print(std::ostream & out, FntExpr<Expr> const & x) { return out<<x; }
 };

 BOOST_PROTO_DEFINE_OPERATORS(is_proto_terminal, FntDomain);

 /* ---------------------------------------------------------------------------------------------------
  *  Does the object have auto_assign
  *  --------------------------------------------------------------------------------------------------- */
#define TRIQS_CLEF_HAS_AUTO_ASSIGN() typedef void __triqs_nvl_has_auto_assign_tag;

 /*
 struct not_defined{};
 template<typename T, typename F> not_defined triqs_nvl_auto_assign(T, F) { return not_defined();};

 template<typename T, typename F> struct has_auto_assign {
  typedef typename boost::unwrap_reference<T>::type T1;
  static T1 & t1();
 // typedef decltype (triqs_nvl_auto_assign( T1 &, F)) R; 
  typedef BOOST_TYPEOF_TPL ( triqs_nvl_auto_assign( t1(),(*(pseudo_default_construct<F>()))) ) R;
  //typedef BOOST_TYPEOF_TPL ( triqs_nvl_auto_assign( *(static_cast<T1*>(NULL)),(*(pseudo_default_construct<F>()))) ) R;
  //typedef BOOST_TYPEOF_TPL ( triqs_nvl_auto_assign((*(pseudo_default_construct<T1>())),(*(pseudo_default_construct<F>()))) ) R;
  typedef typename mpl::not_<boost::is_same<R,not_defined> >::type type;
  static const bool value = type::value;
 };  
  
*/
 template<typename T, typename F, typename Void=void> struct has_auto_assign : mpl::false_{};
 template<typename T, typename F> struct has_auto_assign<T,F,typename T::__triqs_nvl_has_auto_assign_tag> : mpl::true_{};
 template <typename T,typename F> struct has_auto_assign<boost::reference_wrapper<T>,F,void > : has_auto_assign<T,F> {};

#define TRIQS_CLEF_HAS_AUTO_ASSIGN_SUBSCRIPT() typedef void __triqs_nvl_has_auto_assign_subscript_tag;
 template<typename T, typename F, typename Void=void> struct has_auto_assign_subscript : mpl::false_{};
 template<typename T, typename F> struct has_auto_assign_subscript<T,F,typename T::__triqs_nvl_has_auto_assign_subscript_tag>:mpl::true_{};
 template<typename T,typename F> struct has_auto_assign_subscript<boost::reference_wrapper<T>,F,void>:has_auto_assign_subscript<T,F>{};

 /* ---------------------------------------------------------------------------------------------------
  *  Implementation of assignement for expressions 
  * --------------------------------------------------------------------------------------------------- */

 // expression of the form f(z_)(x_,y_) = value
 // x_, y_ MUST be just one place holder and f must be such that is_assignable_leaf<f> = true

 // no compile time check since has_auto_assign could be defined for some F only ....
 //template<typename T> struct is_assignable_leaf : mpl::or_< has_auto_assign_subscript<T>, has_auto_assign<T> > {};
 //template<typename T> struct is_assignable_leaf<boost::reference_wapper<T> > : is_assignable_leaf<T> {};
 template<typename T> struct is_assignable_leaf : mpl::not_< is_lazy <T> > {}; 

 struct AssignableLeafGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<is_assignable_leaf<proto::_value>()> > {}; 

 struct lhs_Grammar;
 struct lhs_call_Grammar : 
  proto::or_<
  proto::function <AssignableLeafGrammar,proto::vararg<PlaceholderLeafGrammar > > 
  ,proto::function <lhs_Grammar,proto::vararg<PlaceholderLeafGrammar > > > {} ;

 struct lhs_subscript_Grammar : 
  proto::or_<
  proto::subscript <AssignableLeafGrammar,PlaceholderLeafGrammar > 
  ,proto::subscript <lhs_Grammar,PlaceholderLeafGrammar  > > {} ;

 struct lhs_Grammar : proto::or_<lhs_call_Grammar,lhs_subscript_Grammar>{};

 template <bool is_call, typename LHS, typename RHS, typename Enable=void> struct assignment_impl { 
  static void invoke(LHS & lhs, RHS const & rhs) { boost::unwrap_ref(lhs) = rhs;} }; 
 template <typename LHS, typename RHS> void assignment_subscript(LHS lhs, RHS rhs){assignment_impl<false,LHS,RHS>::invoke(lhs,rhs);} 
 template <typename LHS, typename RHS> void assignment_call(LHS lhs, RHS rhs){assignment_impl<true,LHS,RHS>::invoke(lhs,rhs);} 

 template <typename LHS, typename RHS>  // objects which have set_from_function
  struct assignment_impl<true,LHS, RHS, typename enable_if< has_auto_assign <LHS,RHS> >::type > { 
   static void invoke(LHS & lhs, RHS const & rhs) {  triqs_nvl_auto_assign(boost::unwrap_ref(lhs), rhs);} 
  }; 

 template <typename LHS, typename RHS>  // objects which have set_from_function
  struct assignment_impl<false,LHS, RHS, typename enable_if< has_auto_assign_subscript <LHS,RHS> >::type > { 
   static void invoke(LHS & lhs, RHS const & rhs) {  triqs_nvl_auto_assign_subscript (boost::unwrap_ref(lhs), rhs);} 
  }; 

 template < typename Expr, typename RHS>  // expressions which are terminal, just undress the terminal and recall
  struct assignment_impl<false,FntExpr<Expr>, RHS, typename enable_if< proto::matches<FntExpr<Expr> , proto::terminal<proto::_> > >::type > 
  { static void invoke(FntExpr<Expr> const & lhs, RHS const & rhs) { triqs::clef::assignment_subscript ( proto::value(lhs), rhs); } };

 template <typename Expr, typename RHS>  // expressions which are terminal, just undress the terminal and recall
  struct assignment_impl<true,FntExpr<Expr>,RHS,typename enable_if< proto::matches<FntExpr<Expr> , proto::terminal<proto::_> > >::type > 
  { static void invoke(FntExpr<Expr> const & lhs, RHS const & rhs) { triqs::clef::assignment_call ( proto::value(lhs), rhs); } };

 template <bool is_call, typename Expr, typename RHS> // expression which are not terminal. just call = again !
  struct assignment_impl<is_call,FntExpr<Expr>,RHS,typename enable_if<mpl::not_<proto::matches<FntExpr<Expr>,proto::terminal<proto::_> > > >::type> 
  { static void invoke(FntExpr<Expr> const & lhs, RHS const & rhs) { lhs = rhs; } }; 

 template <typename Expr, typename T> 
  struct FntExpr_assign_delegation<Expr, T, typename enable_if<proto::matches<FntExpr<Expr> , lhs_Grammar > >::type > {
   typedef void result_type;
   static const bool is_call = proto::matches<FntExpr<Expr> , lhs_call_Grammar >::value;
   static const int n_var = proto::arity_of<FntExpr<Expr> >::value -1;
   static result_type invoke (FntExpr<Expr> expr, T rhs) { invoke_impl(mpl::bool_<is_call>(),mpl::int_<n_var>(), expr,rhs);}
   private : // implementation 
   static result_type invoke_impl (mpl::false_, mpl::int_<1>, FntExpr<Expr> expr, T rhs) {
    triqs::clef::assignment_subscript ( proto::left(expr)  ,make_function(rhs, proto::value(proto::child_c<1>(expr)))); 
   }
#define AUX1(z,p,unused) proto::value(proto::child_c<p+1>(expr) )
#define IMPL(z, NN, unused) \
   static result_type invoke_impl (mpl::true_, mpl::int_<NN>, FntExpr<Expr> expr, T rhs) {\
    triqs::clef::assignment_call ( proto::left(expr)  ,make_function(rhs, BOOST_PP_ENUM(NN,AUX1,nil) ) ); \
   }
   BOOST_PP_REPEAT_FROM_TO(1,BOOST_PP_INC(TRIQS_CLEF_MAXNARGS) , IMPL, nil);
#undef AUX1
#undef IMPL 
  };

 /* ---------------------------------------------------------------------------------------------------
  *   transformation of a couple ( FntExpr, placeholder_list)  into a regular function
  *   Indeed, once we know the order of the placeholder, we can make an ordinary function with positional ordering
  *    
  *    auto make_function( ExpressionType expr, PH1, PH2, PH3, .....);
  * --------------------------------------------------------------------------------------------------- */

 template <typename PH> struct undress_placeholder_expr  {//recover PH from the corresponding proto terminal expression
  static_assert( (proto::matches< PH , PlaceholderLeafGrammar >::value), "Internal Error");
  typedef typename proto::result_of::value< PH >::value_type type; // extract the placeholder type
 };
 template <int N> struct undress_placeholder_expr<placeholder<N> >  { typedef placeholder<N> type; };

 template <typename FE, typename X_vec, int N = mpl::size<X_vec>::value > struct make_function_impl;

#define AUX_TY(z,p,unused)       typedef typename mpl::at_c<X_vec,p>::type X##p;
#define AUX_PAIR(z,p,unused)     bf::make_pair< X##p >(v##p)
#define AUX_PAIR_DEF(z,p,unused) bf::pair< X##p , V##p> 
#define IMPL(z, NN, unused)  \
 template <typename FE, typename X_vec>\
 struct make_function_impl<FE,X_vec,NN> {\
  typedef X_vec ph_vec_exclude;\
  typedef FE expression_type; FE expr; make_function_impl ( FE const & expr_) : expr(expr_) {}\
  BOOST_PP_REPEAT(NN,AUX_TY,nil)\
  template<BOOST_PP_ENUM_PARAMS(NN,typename V)> struct result_type  {\
   typedef bf::map< BOOST_PP_ENUM(NN,AUX_PAIR_DEF,nil) > map_type;\
   typedef typename result_of::eval_on_map<FE, map_type>::type type;\
  };\
  template<BOOST_PP_ENUM_PARAMS(NN,typename V)>\
  typename result_type<BOOST_PP_ENUM_PARAMS(NN, V)>::type\
  operator()  ( BOOST_PP_ENUM_BINARY_PARAMS(NN,V,v)) const {\
   typename result_type<BOOST_PP_ENUM_PARAMS(NN, V)>::map_type m ( BOOST_PP_ENUM(NN,AUX_PAIR,nil)  ); \
   return triqs::clef::eval (expr,BOOST_PP_ENUM(NN,AUX_PAIR,nil) ); \
  }\
  friend std::ostream &operator <<(std::ostream &sout, make_function_impl const & expr) { \
   return sout << " function : ("<< placeholder_list_as_string1<typename make_function_impl::ph_vec_exclude>() <<") --> "<<  expr.expr ; }\
  friend std::ostream & triqs_nvl_formal_print(std::ostream & out, make_function_impl const & x) { return out<<x; }\
 };
 BOOST_PP_REPEAT_FROM_TO(1,BOOST_PP_INC(TRIQS_CLEF_MAXNARGS), IMPL, nil);
#undef AUX_TY
#undef AUX_PAIR_DEF
#undef AUX_PAIR
#undef IMPL 

 /* ---------------------------------------------------------------------------------------------------
  *   transformation of a couple ( FntExpr, placeholder_list)  into a regular function
  *   Indeed, once we know the order of the placeholder, we can make an ordinary function with positional ordering
  * --------------------------------------------------------------------------------------------------- */

 template <typename Expr, typename X_vec>
  struct lazy_function_impl  { 
   lazy_function_impl ( Expr const & expr_) : expr(expr_) {} 

   TRIQS_CLEF_IS_EXPRESSION();   
   typedef Expr expression_type;
   Expr expr;

   struct _l { template <typename T> struct apply : undress_placeholder_expr< T > {}; };
   typedef typename mpl::transform_view <X_vec, _l >::type ph_set_exclude;
   typedef typename mpl::filter_view<typename get_ph_set<Expr>::type,mpl::not_<mpl::contains<ph_set_exclude,mpl::_> > >::type ph_set;

   static const bool evaluated = mpl::empty< ph_set >::value;
   typedef typename  mpl::if_c< evaluated,  make_function_impl <Expr, X_vec >, lazy_function_impl >::type type2;

   template<typename BfMapType> struct call_rtype {
    typedef typename result_of::eval_on_map<Expr,BfMapType>::type R;
    typedef lazy_function_impl <R, X_vec> L;
    typedef typename mpl::if_c< L::evaluated, make_function_impl <R, X_vec >, L>::type type;
   };

   template<typename BfMapType> typename call_rtype<BfMapType>::type  
    eval_on_map (BfMapType const & ph_value_dict) const { return typename call_rtype<BfMapType>::type(expr.eval_on_map(ph_value_dict));}

   friend std::ostream &operator <<(std::ostream &sout, lazy_function_impl const & expr) { 
    return sout << "lazy function : ("<< placeholder_list_as_string1<lazy_function_impl::ph_set_exclude>() <<") --> "<<  expr.expr ; 
   }
   friend std::ostream & triqs_nvl_formal_print(std::ostream & out, lazy_function_impl const & x) { return out<<x; }
  };

 namespace result_of { 
  template<typename Expr, typename X_vec> struct make_function { typedef typename lazy_function_impl<Expr, X_vec >::type2 type;};
 }

#define IMPL(z, N, unused)  \
 template <typename Expr, typename_A(N) > \
 typename result_of::make_function<Expr, mpl::vector<_A_(N)> >::type \
 make_function(Expr expr, _A_(N) ) { return typename lazy_function_impl<Expr,mpl::vector<_A_(N)> >::type2 (expr);}
 BOOST_PP_REPEAT_FROM_TO(1,BOOST_PP_INC(TRIQS_CLEF_MAXNARGS), IMPL, nil);
#undef IMPL 

 /* --------------------------------------------------------------------------------------------------
  *  x_ >> expression  is the same as make_function(expression,x)
  * --------------------------------------------------------------------------------------------------- */

 template <int N, typename E2>
  typename result_of::make_function<typename proto::result_of::as_expr<E2 const,FntDomain >::type, mpl::vector<placeholder<N> > >::type 
  operator >> ( placeholder<N> const & x, E2 const & expr) { return make_function (proto::as_expr<FntDomain>(expr), x); } 

 /* --------------------------------------------------------------------------------------------------
  * lazy (x) make any object lazy
  * --------------------------------------------------------------------------------------------------- */

 template<typename T> 
  typename proto::result_of::as_expr<T const ,FntDomain>::type 
  lazy_copy(T const & x) { return proto::as_expr<FntDomain>(x);}

 template<typename T> 
  typename proto::result_of::as_expr<boost::reference_wrapper<T> const ,FntDomain>::type 
  lazy(T & x) { return proto::as_expr<FntDomain>(boost::ref(x));}

 template<typename T> 
  typename proto::result_of::as_expr<boost::reference_wrapper<T const> const ,FntDomain>::type 
  lazy(T const & x) { return proto::as_expr<FntDomain>(boost::cref(x));}

  template< typename Obj> struct lazy_call {
   Obj const & obj; lazy_call(Obj const & obj_): obj(obj_) {}
#define IMPL(z,N,unused)\
  template< typename_A(N) >\
  typename proto::result_of::make_expr< p_tag::function, FntDomain, typename proto::result_of::as_expr<Obj ,FntDomain>::type,_A_const_ref(N)>::type\
  operator() (_A_const_ref_a(N) ) { return proto::make_expr<p_tag::function, FntDomain>(proto::as_expr<FntDomain>(obj),_a_(N));} 
  BOOST_PP_REPEAT_FROM_TO(1,TRIQS_CLEF_MAXNARGS, IMPL, nil);
#undef IMPL
  };
 
  namespace result_of {

  template< typename Obj, typename Arg > struct subscript_enabled_if_arg_is_lazy : 
   enable_if< is_lazy<Arg> , 
   typename proto::result_of::make_expr< boost::proto::tag::subscript, FntDomain, Obj, Arg const &>::type > {};

  template< typename Obj, BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT (TRIQS_CLEF_MAXNARGS,typename Arg, void ) > 
   struct call_enabled_if_one_arg_is_lazy;

#define IMPL(z,N,unused)\
  template< typename Obj, typename_A(N) >\
  struct call_enabled_if_one_arg_is_lazy<Obj,_A_void(N) > :\
  enable_if< one_is_lazy<_A_(N)>,\
  typename proto::result_of::make_expr< p_tag::function, FntDomain, typename proto::result_of::as_expr<Obj ,FntDomain>::type,_A_const_ref(N)>::type\
  >{};
  BOOST_PP_REPEAT_FROM_TO(1,TRIQS_CLEF_MAXNARGS, IMPL, nil);
#undef IMPL

  // not ok on 4.6
  /*template< typename Obj, typename... Args >
  struct call_enabled_if_one_arg_is_lazy_v: enable_if< one_is_lazy<Args...>,
  typename proto::result_of::make_expr< p_tag::function, FntDomain, typename proto::result_of::as_expr<Obj ,FntDomain>::type,Args const & ...>::type
  >{};
 */
 
 }
}} //  namespace triqs::clef
/* *****************************************************************************************
 *    MACROS 
 * *****************************************************************************************/
/* --------------------------------------------------------------------------------------------------
 *  The macro to add the () operator on expressions for any callable object
 *  TRIQS_CLEF_ADD_LAZY_CALL_WITH_COPY(Number_of_arguments, thetype_being_defined) : put in the definition of the object to produce the () operator catching the placeholder
 *  TRIQS_CLEF_ADD_LAZY_CALL_REF(Number_of_arguments, thetype_being_defined) : put in the definition of the object to produce the () operator catching the placeholder
 *
 *  TRIQS_CLEF_ADD_LAZY_SUBSCRIPT(thetype_being_defined) : put in the definition of the object to produce the [] operator catching the placeholder
 *  TRIQS_CLEF_ADD_LAZY_SUBSCRIPT_REF(thetype_being_defined) : put in the definition of the object to produce the [] operator catching the placeholder
 *
 *  The first version forces a copy of the object when called, the second takes a reference.
 * --------------------------------------------------------------------------------------------------- */

#define TRIQS_CLEF_ADD_LAZY_SUBSCRIPT_WITH_COPY(MYTYPE)\
 template< typename Arg > \
const typename triqs::clef::result_of::subscript_enabled_if_arg_is_lazy<MYTYPE, Arg>::type  \
operator[]( Arg const & arg ) const { return triqs::clef::lazy_copy(*this)[arg];}

#define TRIQS_CLEF_ADD_LAZY_SUBSCRIPT_REF(MYTYPE)\
 template< typename Arg > \
const typename triqs::clef::result_of::subscript_enabled_if_arg_is_lazy<boost::reference_wrapper<MYTYPE const>, Arg>::type  \
operator[]( Arg const & arg ) const { return triqs::clef::lazy(*this)[arg];}\
template< typename Arg > \
const typename triqs::clef::result_of::subscript_enabled_if_arg_is_lazy<boost::reference_wrapper<MYTYPE>, Arg>::type  \
operator[]( Arg const & arg ) { return triqs::clef::lazy(*this)[arg];}

#define TRIQS_CLEF_ADD_LAZY_CALL_WITH_COPY(P,MYTYPE) \
 template< BOOST_PP_ENUM_PARAMS(P,typename Arg) >\
const typename triqs::clef::result_of::call_enabled_if_one_arg_is_lazy<MYTYPE, BOOST_PP_ENUM_PARAMS(P,Arg)>::type  \
operator()( BOOST_PP_ENUM_BINARY_PARAMS(P,Arg,const & arg) ) const { return triqs::clef::lazy_copy(*this)(BOOST_PP_ENUM_PARAMS(P,arg));}

#define TRIQS_CLEF_ADD_LAZY_CALL_WITH_VIEW(P,VIEWTYPE) \
 template< BOOST_PP_ENUM_PARAMS(P,typename Arg) >\
const typename triqs::clef::result_of::call_enabled_if_one_arg_is_lazy<VIEWTYPE, BOOST_PP_ENUM_PARAMS(P,Arg)>::type  \
operator()( BOOST_PP_ENUM_BINARY_PARAMS(P,Arg,const & arg) ) const { return triqs::clef::lazy_copy(VIEWTYPE(*this))(BOOST_PP_ENUM_PARAMS(P,arg));}

#define TRIQS_CLEF_ADD_LAZY_CALL_REF(P,MYTYPE) \
 template< BOOST_PP_ENUM_PARAMS(P,typename Arg) >\
const typename triqs::clef::result_of::call_enabled_if_one_arg_is_lazy<boost::reference_wrapper<MYTYPE const>, BOOST_PP_ENUM_PARAMS(P,Arg)>::type  \
operator()( BOOST_PP_ENUM_BINARY_PARAMS(P,Arg,const & arg) ) const { return triqs::clef::lazy(*this)(BOOST_PP_ENUM_PARAMS(P,arg));}\
template< BOOST_PP_ENUM_PARAMS(P,typename Arg) >\
const typename triqs::clef::result_of::call_enabled_if_one_arg_is_lazy<boost::reference_wrapper<MYTYPE>, BOOST_PP_ENUM_PARAMS(P,Arg)>::type  \
operator()( BOOST_PP_ENUM_BINARY_PARAMS(P,Arg,const & arg) ) { return triqs::clef::lazy(*this)(BOOST_PP_ENUM_PARAMS(P,arg));}

/* --------------------------------------------------------------------------------------------------
 *  The macro to make any function lazy
 *  TRIQS_CLEF_MAKE_FNT_LAZY_WRAP (Arity,FunctionName ) : creates a new function in the triqs::lazy namespace 
 *  taking expressions (at least one argument has to be an expression) and returning the proto object.
 *  The lookup happens by ADL, so IT MUST BE USED IN THE triqs::lazy namespace
 * --------------------------------------------------------------------------------------------------- */

#define TRIQS_CLEF_MAKE_FNT_LAZY_aux_aux(z,P,unused) A##P()
#define TRIQS_CLEF_MAKE_FNT_LAZY_aux(ARITY,name)\
 struct name##_lazy_impl {\
  template<BOOST_PP_ENUM_PARAMS(ARITY,typename A) > struct _res { typedef BOOST_TYPEOF_TPL(name(BOOST_PP_ENUM(ARITY,TRIQS_CLEF_MAKE_FNT_LAZY_aux_aux,nil))) type;};\
  template< BOOST_PP_ENUM_PARAMS(ARITY,typename A) > typename _res<BOOST_PP_ENUM_PARAMS(ARITY,A) >::type operator() (BOOST_PP_ENUM_BINARY_PARAMS(ARITY,A,const & a) ) const \
  { return name(BOOST_PP_ENUM_PARAMS(ARITY,a));}\
 };

#define TRIQS_CLEF_MAKE_FNT_LAZY_NOWRAP(P,NAME, TY, FNT )\
 template< BOOST_PP_ENUM_PARAMS(P,typename Arg) >\
const typename triqs::clef::result_of::call_enabled_if_one_arg_is_lazy<TY, BOOST_PP_ENUM_PARAMS(P,Arg)>::type  \
NAME( BOOST_PP_ENUM_BINARY_PARAMS(P,Arg,const & arg) ){ return triqs::clef::lazy_copy(FNT)(BOOST_PP_ENUM_PARAMS(P,arg));}

#define TRIQS_CLEF_MAKE_FNT_LAZY(Arity,name) TRIQS_CLEF_MAKE_FNT_LAZY_aux(Arity,name) TRIQS_CLEF_MAKE_FNT_LAZY_NOWRAP(Arity,name, name##_lazy_impl,name##_lazy_impl()) 

#undef ALLOWED_BINARY_OP
#undef OP_OP
#undef OP_NAME
#undef typename_A
#undef _A_const_ref
#undef _A_const_ref_a
#undef _A_const_ref_a_tr
#undef _A_ref
#undef _A_
#undef _a_
#undef _A_void
#endif 
