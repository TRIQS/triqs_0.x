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
#ifndef TRIQS_ARRAYS_EXPRESSION_ARRAY_ALGEBRA_H
#define TRIQS_ARRAYS_EXPRESSION_ARRAY_ALGEBRA_H
#include <triqs/utility/proto/tools.hpp>
#include "../array.hpp"
namespace triqs { namespace arrays {

 namespace tup = triqs::utility::proto; namespace mpl = boost::mpl; namespace proto = boost::proto;

 template<typename Expr> struct array_expr;

 namespace detail_expr_array { 

  using boost::enable_if; using boost::enable_if_c; using proto::_left; using proto::_right; 
  using proto::_value; namespace p_tag= proto::tag; using boost::remove_reference;

  struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<tup::is_in_ZRC<proto::_value>()> > {}; 
  struct BasicArrayTypeGrammar :  proto::and_< proto::terminal<proto::_>, proto::if_<ImmutableCuboidArray<proto::_value>()> > {}; 

  struct ArrayGrammar : 
   proto::or_<
   ScalarGrammar // scalar in injected in the algebra 
   , BasicArrayTypeGrammar
   , proto::plus      <ArrayGrammar,ArrayGrammar>
   , proto::minus     <ArrayGrammar,ArrayGrammar>
   , proto::multiplies<ArrayGrammar,ArrayGrammar>
   , proto::divides   <ArrayGrammar,ArrayGrammar>
   , proto::negate    <ArrayGrammar >
   > {};

 struct eval_scalar { // a transform that evaluates a scalar as an array with the scalar at all position
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename S, typename KeyType> struct result<This(S,KeyType)> { typedef typename remove_reference<S>::type type;};
   template<typename S, typename KeyType> 
    typename result<eval_scalar(S,KeyType)>::type operator ()(S const & s, KeyType const & key) const {return s;}
  };

  struct eval_mv { // a transform that evaluates an array on a Key 
   BOOST_PROTO_CALLABLE();

   template<typename Sig> struct result;
   template<typename This, typename M, typename KeyType> struct result<This(M,KeyType)> { typedef typename remove_reference<M>::type::value_type type;};
   template<typename M, typename KeyType> 
    typename result<eval_mv(M,KeyType)>::type operator ()(M const & m, KeyType const & key) const { return m[key]; }
  };

  struct eval_t;
  struct eval_t_cases { template <typename TAG> struct case_: proto::not_<proto::_>{}; };
  template<> struct eval_t_cases::case_<proto::tag::terminal> :  
   proto::or_<
   proto::when<ScalarGrammar, eval_scalar(proto::_value,proto::_state)>
   ,proto::when<BasicArrayTypeGrammar, eval_mv (proto::_value, proto::_state) >
   > {};
  template<> struct eval_t_cases::case_<proto::tag::plus>  : proto::when< proto::plus <eval_t,eval_t>,  proto::_default<eval_t> > {};
  template<> struct eval_t_cases::case_<proto::tag::minus> : proto::when< proto::minus <eval_t,eval_t>, proto::_default<eval_t> > {};
  template<> struct eval_t_cases::case_<proto::tag::multiplies> : proto::when< proto::minus <eval_t,eval_t>, proto::_default<eval_t> > {};
  template<> struct eval_t_cases::case_<proto::tag::divides> : proto::when< proto::minus <eval_t,eval_t>, proto::_default<eval_t> > {};
  template<> struct eval_t_cases::case_<proto::tag::negate> : proto::when< proto::negate<eval_t>, tup::negate_t(eval_t(_left))> {};
  struct eval_t : proto::switch_<eval_t_cases> {};

  // -----------  computation of the domain -------------------

  struct no_domain {static const int rank = 0; typedef mini_vector<size_t,0> index_value_type;};

  struct get_domain {
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename X> struct result<This(X)> { typedef typename remove_reference<X>::type::domain_type type;};
   template<typename X> typename X::domain_type operator ()(X const & x) const { return x.domain();}
  };

  struct combine_domain {
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename D> struct result<This(D,D)>          {typedef D type;};
   template<typename This, typename D> struct result<This(D,no_domain)>  {typedef D type;};
   template<typename This, typename D> struct result<This(no_domain,D)>  {typedef D type;};
   template<typename D> D operator ()(no_domain const & d1, D const & d2) const { return d2;}
   template<typename D> D operator ()(D const & d1, no_domain const & d2) const { return d1;}
   template<typename D> D operator ()(D const & d1, D const & d2) const { 
    if (d1.lengths() != d2.lengths()) TRIQS_RUNTIME_ERROR << "Domain size mismatch : "<< d1.lengths()<<" vs" <<d2.lengths();
    return d1;
   } 
  };

  struct dom_t : 
   proto::or_<
   proto::when< ScalarGrammar, no_domain() >
   ,proto::when< BasicArrayTypeGrammar, get_domain(proto::_value) >
   ,proto::when< proto::binary_expr <proto::_,dom_t,dom_t>,  combine_domain (dom_t(proto::_left), dom_t( proto::_right)) >
   ,proto::when< proto::unary_expr<proto::_,dom_t >, dom_t(proto::_left) >
   > {};

 /* ---------------------------------------------------------------------------------------------------
   * Define the main expression template ArrayExpr, the domain and Grammar (implemented below)
   * NB : it modifies the PROTO Domain to make *COPIES* of all objects.
   * We escape the copy of array by specializing the template below (end of file)
   * ---> to be rediscussed.
   * cf http://www.boost.org/doc/libs/1_49_0/doc/html/proto/users_guide.html#boost_proto.users_guide.front_end.customizing_expressions_in_your_domain.per_domain_as_child
   --------------------------------------------------------------------------------------------------- */
  struct ArrayDomain : proto::domain<proto::generator<array_expr>, ArrayGrammar> {
   template< typename T > struct as_child : proto_base_domain::as_expr< T > {};
  };

  //For arrays special treatment : the regular classes are replaced by the corresponding const view
  template< typename T, int N ,typename Opt> struct ArrayDomain::as_child< array<T,N,Opt> > : 
   ArrayDomain::proto_base_domain::template as_expr< const array_view<T,N,Opt> >{};

  template< typename T, int N, typename Opt> struct ArrayDomain::as_child< const array<T,N,Opt> > : 
   ArrayDomain::proto_base_domain::template as_expr< const array_view<T,N,Opt> >{};

 } // namespace detail_expr_array

 //   Expression for arrays
 template<typename Expr> struct array_expr : 
  TRIQS_MODEL_CONCEPT(ImmutableCuboidArray),       proto::extends<Expr, array_expr<Expr>, detail_expr_array::ArrayDomain> { 
  array_expr( Expr const & expr = Expr() ) : proto::extends<Expr, array_expr<Expr>, detail_expr_array::ArrayDomain> ( expr ) {} 

  typedef size_t index_type;
  typedef typename boost::remove_reference<typename boost::result_of<detail_expr_array::dom_t(array_expr) >::type>::type domain_type;
  static const int rank = domain_type::rank; 
  typedef typename domain_type::index_value_type key_type; 
  typedef typename boost::remove_reference<typename boost::result_of<detail_expr_array::eval_t(array_expr,key_type) >::type>::type value_type;

  domain_type domain() const { return detail_expr_array::dom_t()(*this); }
  mini_vector<size_t,2> shape() const { return this->domain().lengths();} 
  value_type operator[] (key_type const &key) const { return detail_expr_array::eval_t()(*this, key); }

  friend std::ostream &operator <<(std::ostream &sout, array_expr<Expr> const &expr){return proto::eval(expr,tup::algebra_print_ctx (sout));}
 };

 template<typename Expr> struct ImmutableCuboidArray<array_expr<Expr> > : 
  mpl::not_<boost::proto::matches<Expr,detail_expr_array::ScalarGrammar> >{}; // every expression but a single scalar

 BOOST_PROTO_DEFINE_OPERATORS(ImmutableCuboidArray, detail_expr_array::ArrayDomain);

 template<typename Expr > array_view <typename Expr::value_type, Expr::domain_type::rank>
  make_array( Expr const & e) { return array<typename Expr::value_type, Expr::domain_type::rank>(e);}

}}//namespace triqs::arrays

#endif
