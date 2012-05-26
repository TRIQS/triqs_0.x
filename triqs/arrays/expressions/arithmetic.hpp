
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

#ifndef TRIQS_ARRAYS_EXPRESSION_ARITH_H
#define TRIQS_ARRAYS_EXPRESSION_ARITH_H

#include <boost/mpl/int.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/proto/core.hpp>
#include <boost/proto/context.hpp>
#include <boost/proto/transform/arg.hpp>
#include <boost/proto/transform.hpp>
#include "../../utility/typeid_name.hpp"
#include "../array.hpp"
#include "../matrix.hpp"

namespace mpl = boost::mpl;
namespace proto = boost::proto;
namespace triqs { namespace arrays { namespace expressions { namespace arithmetic { 

 struct ScalarGrammar : 
  proto::or_<
  proto::terminal< int >
  , proto::terminal< long >
  , proto::terminal< float >
  , proto::terminal< double >
  , proto::terminal< std::complex<float> >
  , proto::terminal< std::complex<double> >
  >
 {};

 //template<typename T> struct IsArray : Tag::check<Tag::expression_terminal,T> {}; 

 struct BasicArrayTypeGrammar :  
  proto::and_< 
  proto::terminal<proto::_>, 
  proto::if_<ImmutableCuboidArray<proto::_value>()> 
  >
 {}; 

 struct ArrayGrammar : 
  proto::or_<
  ScalarGrammar
  , BasicArrayTypeGrammar
  , proto::plus< ArrayGrammar, ArrayGrammar >
  , proto::minus< ArrayGrammar, ArrayGrammar >
  , proto::multiplies< ScalarGrammar, ArrayGrammar >
  , proto::multiplies< ArrayGrammar, ScalarGrammar >
  , proto::divides< ArrayGrammar, ScalarGrammar >
  >
 {};

 // given 2 domains and an op, returns the domain of the result of op
 template<typename ProtoTag, typename D1,typename D2> struct _dom_ops;

 template<typename D1, typename D2> struct _dom_ops<proto::tag::plus,D1,D2> { 
  static_assert ( (boost::is_same<D1,D2>::value), "Error : adding/substracting 2 arrays with different ranks !");
  typedef D1 domain_type;
  static domain_type invoke(D1 const & a, D2 const & b) { 
   TRIQS_ARRAYS_DEBUG_CHECK((a==b), "Error : adding/substracting 2 arrays with different domains"<<a<<" "<<b);
   return a;
  }
 }; 

 template<typename D1, typename D2> struct _dom_ops<proto::tag::multiplies,D1,D2> {
  typedef typename triqs::arrays::indexmaps::result_of::tensor_product<D1,D2>::type domain_type;
  static domain_type invoke(D1 const & d1, D2 const & d2) {return tensor_product(d1,d2);} 
 };

 template<typename D1, typename D2> struct _dom_ops<proto::tag::minus,D1,D2>:  _dom_ops<proto::tag::plus,D1,D2> {};
 template<typename D1, typename D2> struct _dom_ops<proto::tag::divides,D1,D2>:_dom_ops<proto::tag::multiplies,D1,D2> {};

 template<typename ProtoTag, typename T1, typename T2> struct TypeDomain { 
  typedef typename T1::value_type V1;
  typedef typename T2::value_type V2;
  typedef BOOST_TYPEOF_TPL( V1() + V2()) value_type;
  typedef typename _dom_ops< ProtoTag, typename T1::domain_type, typename T2::domain_type>::domain_type domain_type;
  T1 const & a; T2 const & b;
  TypeDomain (T1 const & a_, T2 const & b_):a(a_),b(b_) {} 
  domain_type domain() const { return _dom_ops<ProtoTag,typename T1::domain_type, typename T2::domain_type>::invoke(a.domain(),b.domain()); }
 }; 

 template<typename T> struct wrap_scalar { 
  typedef T value_type;
  typedef indexmaps::cuboid_domain<0> domain_type;
  wrap_scalar( T const &) {}
  domain_type domain() const { return indexmaps::cuboid_domain<0>();}
 };

 // the main proto transform used here...
 using proto::_left; using proto::_right;
 struct Domain : 
  proto::or_<
  proto::when< ScalarGrammar, wrap_scalar<proto::_value>(proto::_value) >
  ,proto::when< proto::terminal< proto::_ > , proto::_value >
  ,proto::when< proto::plus <Domain, Domain >, TypeDomain<proto::tag::plus, _left, _right > (_left,_right) >
  ,proto::when< proto::minus <Domain, Domain >, TypeDomain<proto::tag::minus, _left, _right > (_left,_right)  >
  ,proto::when< proto::multiplies<Domain, Domain >, TypeDomain< proto::tag::multiplies, _left, _right > (_left,_right)  >
  ,proto::when< proto::divides<Domain, Domain >, TypeDomain< proto::tag::divides, _left, _right > (_left,_right)  >
  ,proto::when< proto::unary_expr<proto::_, Domain >, Domain(proto::_child0)  >
  >
 {};

 template<typename Expr> struct ArrayExpr;
 struct ArrayDomain : proto::domain<proto::generator<ArrayExpr>, ArrayGrammar> {};

 /* -------------------------------------------
  *   Evaluation context
  * ------------------------------------------ */

 template<typename KeyType, typename ReturnType>
  struct ArrayEvalCtx : proto::callable_context< ArrayEvalCtx<KeyType,ReturnType> const > {
   typedef ReturnType result_type;
   KeyType const & key;
   ArrayEvalCtx(KeyType const & key_) : key(key_) {}
   template<typename T>// overrule just the terminals which have array interface.
    typename boost::enable_if< ImmutableArray<T>, result_type >::type 
    operator ()(proto::tag::terminal, T const & t) const { return t[key]; }
  };

 /* -------------------------------------------
  *   Print context
  * ------------------------------------------ */
 template <typename T> std::ostream & formal_print(std::ostream & out, T const & x) { return out<<x;}

 // Formal print of A, no data.
 struct ArrayPrintCtx : proto::callable_context< ArrayPrintCtx const > {
  typedef std::ostream &result_type;
  result_type out;
  ArrayPrintCtx(std::ostream & out_):out(out_) {}
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

 /* -------------------------------------------
  *   Expression
  * ------------------------------------------ */
 template<typename Expr>
  struct ArrayExpr : TRIQS_MODEL_CONCEPT(ImmutableCuboidArray), proto::extends<Expr, ArrayExpr<Expr>, ArrayDomain> {
   typedef proto::extends<Expr, ArrayExpr<Expr>, ArrayDomain> base_type;
   ArrayExpr( Expr const & expr = Expr() ) : base_type( expr ) {}

   // implements HasImmutableArrayInterface 
   typedef typename boost::result_of<Domain(Expr) >::type _T;
   typedef typename _T::value_type value_type;
   typedef typename _T::domain_type domain_type;
   typedef typename domain_type::index_value_type key_type;
   typedef size_t index_type;
   static const int rank = _T::domain_type::rank;

   domain_type domain() const { return Domain()(*this).domain(); }
   value_type operator[] (key_type const &key) const { return proto::eval(*this, ArrayEvalCtx<key_type ,value_type> (key)); }

   // todo : generalize with boost_PP
   value_type operator()(index_type const & i1) const { 
    static_assert( (rank==1), "rank must be 1");
    return (*this)[mini_vector<size_t,1>(i1)];
   }

   value_type operator()(index_type const & i1, index_type const & i2) const { 
    static_assert( (rank==2), "rank must be 2");
    return (*this)[mini_vector<size_t,2>(i1,i2)];
   }

   std::ostream & print_structure (std::ostream &sout, bool recursive=false, int sh=0) const { 
    sout<<"Expression  "<<*this <<std::endl;//Tools::typeid_name(*this)<<std::endl;
    sout<<"            : return_type : "<<triqs::utility::typeid_name((value_type*)0)<<std::endl;
    sout<<"            : domain_type : "<<triqs::utility::typeid_name((domain_type*)0)<<std::endl;
    return sout;
   }

   friend std::ostream &operator <<(std::ostream &sout, ArrayExpr<Expr> const &expr) { return proto::eval(expr, ArrayPrintCtx (sout)); }

  };
}}

// This makes array and array_view proto terminals
BOOST_PROTO_DEFINE_OPERATORS(ImmutableCuboidArray, expressions::arithmetic::ArrayDomain);

// C++0x only ...  //template<typename Expr, typename ReturnType= array_view <typename Expr::value_type, Expr::domain_type::rank> >
// ReturnType eval( Expr const & e) { return typename ReturnType::non_view_type(e);}

template<typename Expr >
array_view <typename Expr::value_type, Expr::domain_type::rank>
eval( Expr const & e) { return array<typename Expr::value_type, Expr::domain_type::rank>(e);}

template<typename Expr >
matrix_view <typename Expr::value_type>
eval_as_matrix( Expr const & e) { return matrix<typename Expr::value_type>(e);}

}}//namespace triqs::arrays 

#endif

