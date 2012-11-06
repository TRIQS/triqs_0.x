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

#include "./tools.hpp"

namespace triqs { namespace utility { namespace proto { 

 namespace mpl = boost::mpl; namespace proto = boost::proto; namespace p_tag= proto::tag;

 // technical : evaluation of arithmetic operators
 template<typename TAG, typename A, typename B> struct _binary_ops_;
 template<typename TAG, typename A> struct _unary_ops_;

 // the list of allowed operators in the grammar
#define BINARY_OP_LIST ((plus,+))((minus,-))((multiplies,*))((divides,/))
#define UNARY_OP_LIST  ((negate,-))

#define OP_NAME(elem) BOOST_PP_TUPLE_ELEM(2,0,elem)
#define OP_OP(elem)   BOOST_PP_TUPLE_ELEM(2,1,elem)

 // a trick to avoid putting T() in the typeof type deduction ! This code is NEVER used 
 template<class T> typename boost::unwrap_reference<T>::type pdc() { 
  typename boost::unwrap_reference<T>::type * x= NULL; assert(0); return *x; 
 }

#define AUX(r, data, elem) \
 template<typename A, typename B> struct _binary_ops_<p_tag::OP_NAME(elem),A,B> {\
  typedef BOOST_TYPEOF_TPL( pdc<A>() OP_OP(elem) pdc<B>()) result_type;\
  static result_type invoke( A const &a , B const & b) { return a OP_OP(elem) b;} };
 BOOST_PP_SEQ_FOR_EACH(AUX, nil , BINARY_OP_LIST); 
#undef AUX

#define AUX(r, data, elem) \
 template<typename A> struct _unary_ops_<p_tag::OP_NAME(elem),A> {\
  typedef BOOST_TYPEOF_TPL( OP_OP(elem) pdc<A>()) result_type;\
  static result_type invoke( A const &a ) { return OP_OP(elem) a;} };
 BOOST_PP_SEQ_FOR_EACH(AUX, nil , UNARY_OP_LIST); 
#undef AUX

#undef OP_NAME
#undef OP_OP
#undef BINARY_OP_LIST
#undef UNARY_OP_LIST

 template <typename T, typename A0> struct call_result_type { 
  typedef BOOST_TYPEOF_TPL (pdc<T>() (pdc<A0>())) type;
 };

 namespace algebra { 

  template< typename OpsCompound, template<typename T> class is_element, template<typename T> class is_scalar = is_in_ZRC > 
   struct grammar_generator {

    struct LeafGrammar   : proto::and_< proto::terminal<proto::_>, proto::if_<is_element<proto::_value>()> > {}; 
    struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<is_scalar<proto::_value>()> > {}; 

//#define SWITCH_GRAMMAR
#ifndef SWITCH_GRAMMAR

    struct Grammar : 
     proto::or_<
     proto::when< ScalarGrammar,                       typename OpsCompound::template scalar<proto::_value>(proto::_value) >
     ,proto::when< LeafGrammar,                        proto::_value >
     ,proto::when< proto::plus <Grammar,Grammar>,      typename OpsCompound::template plus<proto::_left,proto::_right > (proto::_left,proto::_right) >
     ,proto::when< proto::minus <Grammar,Grammar>,     typename OpsCompound::template minus<proto::_left,proto::_right > (proto::_left,proto::_right)  >
     ,proto::when< proto::multiplies<Grammar,Grammar>, typename OpsCompound::template multiplies<proto::_left,proto::_right > (proto::_left,proto::_right)>
     ,proto::when< proto::divides<Grammar,Grammar>,    typename OpsCompound::template divides<proto::_left,proto::_right > (proto::_left,proto::_right)>
     ,proto::when< proto::negate<Grammar >,            typename OpsCompound::template negate <proto::_left >(proto::_left) >
     > {};

#else
    struct Grammar;
    struct GrammarCases { template <typename TAG> struct case_: proto::not_<proto::_>{}; };
    struct Grammar : proto::switch_<GrammarCases> {};
//#endif
//   #ifdef SWITCH_GRAMMAR

 // template< typename OpsCompound, template<typename T> class is_element, template<typename T> class is_scalar >   
   //template<> struct grammar_generator<OpsCompound,is_element,is_scalar>::GrammarCases::case_<proto::tag::terminal> :  
   template<> struct GrammarCases::case_<proto::tag::terminal> :  
   proto::or_<
   proto::when< ScalarGrammar,                       typename OpsCompound::template scalar<proto::_value>(proto::_value) >
   ,proto::when< LeafGrammar,                        proto::_value >
   > {};

   template<> struct GrammarCases::case_<proto::tag::plus> :  
   proto::when< proto::plus <Grammar,Grammar>,      typename OpsCompound::template plus<proto::_left,proto::_right > (proto::_left,proto::_right) >
   {};

   template<> struct GrammarCases::case_<proto::tag::minus> :  
   proto::when< proto::minus <Grammar,Grammar>,     typename OpsCompound::template minus<proto::_left,proto::_right > (proto::_left,proto::_right)  >
   {};

   template<> struct GrammarCases::case_<proto::tag::multiplies> :  
   proto::when< proto::multiplies<Grammar,Grammar>, typename OpsCompound::template multiplies<proto::_left,proto::_right > (proto::_left,proto::_right)>
   {};

   template<> struct GrammarCases::case_<proto::tag::divides> :  
   proto::when< proto::divides<Grammar,Grammar>,    typename OpsCompound::template divides<proto::_left,proto::_right > (proto::_left,proto::_right)>
   {};

   template<> struct GrammarCases::case_<proto::tag::negate> :  
   proto::when< proto::negate<Grammar >,            typename OpsCompound::template negate <proto::_left >(proto::_left) >
   {};

#endif
 typedef Grammar type;

   };


  template< typename OpsCompound, template<typename T> class is_element, template<typename T> class is_scalar = is_in_ZRC > 
   struct grammar_generator1 {

    struct LeafGrammar   : proto::and_< proto::terminal<proto::_>, proto::if_<is_element<proto::_value>()> > {}; 
    struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<is_scalar<proto::_value>()> > {}; 
    struct Grammar : 
     proto::or_<
     proto::when< ScalarGrammar,                       typename OpsCompound::template scalar<proto::_value>(proto::_value) >
     ,proto::when< LeafGrammar,                        proto::_value >
     ,proto::when< proto::plus <Grammar,Grammar>,      typename OpsCompound::template binary_node<p_tag::plus,proto::_left,proto::_right > (proto::_left,proto::_right) >
     ,proto::when< proto::minus <Grammar,Grammar>,     typename OpsCompound::template binary_node<p_tag::minus,proto::_left,proto::_right > (proto::_left,proto::_right)  >
     ,proto::when< proto::multiplies<Grammar,Grammar>, typename OpsCompound::template binary_node<p_tag::multiplies,proto::_left,proto::_right > (proto::_left,proto::_right)>
     ,proto::when< proto::divides<Grammar,Grammar>,    typename OpsCompound::template binary_node<p_tag::divides,proto::_left,proto::_right > (proto::_left,proto::_right)>
     ,proto::when< proto::negate<Grammar >,            typename OpsCompound::template negate <proto::_left >(proto::_left) >
     > {};

    typedef Grammar type;
   };

  /* -------------------------------------------
   *  Structure of algebra for algebra valued functions
   * ------------------------------------------ */
  struct algebra_function_desc { 

   template<typename S> struct scalar {
    S s; scalar( S const & x) : s(x) {}
    template <typename T> struct call_rtype { typedef S type; };
    template<typename T> typename call_rtype<T>::type operator() (T) const {return s;}
   };

   template<typename ProtoTag, typename L, typename R> struct binary_node  { 
    L const & l; R const & r; binary_node (L const & l_, R const & r_):l(l_),r(r_) {}
    template <typename T> struct call_rtype {
     typedef _binary_ops_<ProtoTag, typename call_result_type<L,T>::type , typename call_result_type<R,T>::type  > ops_type;
     typedef typename ops_type::result_type type;
    };
    template<typename T> typename call_rtype<T>::type operator() (T const & arg) const {return call_rtype<T>::ops_type::invoke(l(arg),r(arg));}
   }; 

   template<typename L> struct negate  { 
    L const & l; negate (L const & l_):l(l_) {} 
    template <typename T> struct call_rtype {
     typedef _unary_ops_<p_tag::negate, typename call_result_type<L,T>::type > ops_type;
     typedef typename ops_type::result_type type;
    };
    template<typename T> typename call_rtype<T>::type operator() (T const & arg) const {return call_rtype<T>::ops_type::invoke(l(arg));}
   }; 
  };

  /* -------------------------------------------
   *   Print context
   * ------------------------------------------ */

  struct print_ctx : proto::callable_context< print_ctx const > {
   typedef std::ostream &result_type;
   result_type out;
   print_ctx(std::ostream & out_):out(out_) {}
   template <typename T>
    typename boost::enable_if<is_in_ZRC<T>, result_type>::type 
    operator ()(proto::tag::terminal, const T & A) const { return out<<A; }
   template <typename T>
    typename boost::disable_if<is_in_ZRC<T>, result_type>::type 
    operator ()(proto::tag::terminal, const T & A) const { return out<< triqs::utility::typeid_name(A); }
   template<typename L, typename R> result_type operator ()(proto::tag::plus, L const &l, R const &r) const { return out << '(' << l << " + " << r << ')'; }
   template<typename L, typename R> result_type operator ()(proto::tag::minus, L const &l, R const &r) const { return out << '(' << l << " - " << r << ')'; }
   template<typename L, typename R> result_type operator ()(proto::tag::multiplies, L const &l, R const &r) const { return out << l << " * " << r; }
   template<typename L, typename R> result_type operator ()(proto::tag::divides, L const &l, R const &r) const { return out << l << " / " << r; }
  };

 } //namespace algebra



 template<typename Grammar> struct domain_and_expression_generator { 
  template <typename Expr> struct The_Expr;  // the expression
  typedef domain<Grammar,The_Expr,true>  expr_domain; // the domain 
  template<typename Expr> struct The_Expr : boost::proto::extends<Expr, The_Expr<Expr>, expr_domain>{ // impl the expression
   typedef boost::proto::extends<Expr, The_Expr<Expr>, expr_domain> basetype;
   The_Expr( Expr const & expr = Expr() ) : basetype ( expr ) {}
   typedef typename boost::result_of<Grammar(Expr) >::type _G;
   template<typename T> typename call_result_type<_G,T>::type operator() (T const & x) const { return Grammar()(*this)(x); }
   friend std::ostream &operator <<(std::ostream &sout, The_Expr<Expr> const &expr) { return boost::proto::eval(expr, triqs::utility::proto::algebra::print_ctx (sout)); }
  };

 };

#define TRIQS_PROTO_DEFINE_ALGEBRA_VALUED_FNT_ALG(type_identification_trait,scalar_identification_trait)\
 BOOST_PROTO_DEFINE_OPERATORS(type_identification_trait, (triqs::utility::proto::domain_and_expression_generator<\
    triqs::utility::proto::algebra::grammar_generator<triqs::utility::proto::algebra::algebra_function_desc,type_identification_trait, scalar_identification_trait >::type\
    >::expr_domain ));

#define TRIQS_PROTO_DEFINE_ALGEBRA_VALUED_FNT_ALG_WITH_DESC(type_identification_trait,scalar_identification_trait, ALG_DESC)\
 typedef triqs::utility::proto::domain_and_expression_generator<\
 triqs::utility::proto::algebra::grammar_generator<\
 ALG_DESC,type_identification_trait, scalar_identification_trait\
 >::type\
 > _D##type_identification_trait;\
 BOOST_PROTO_DEFINE_OPERATORS(type_identification_trait, _D##type_identification_trait::mydomain );

#define TRIQS_PROTO_DEFINE_ALGEBRA_VALUED_FNT_ALG2(type_identification_trait,scalar_identification_trait)\
 TRIQS_PROTO_DEFINE_ALGEBRA_VALUED_FNT_ALG_WITH_DESC(type_identification_trait,scalar_identification_trait, triqs::utility::proto::algebra::algebra_function_desc)

}}}

#endif

