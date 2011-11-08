
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
#include "./algebra_common.hpp"
#include "../array.hpp"
namespace triqs { namespace arrays { namespace expressions { namespace array_algebra { 

 template<typename T> struct IsArray : Tag::check<Tag::array_algebra_expression_terminal,T> {}; 
 struct BasicArrayTypeGrammar :  proto::and_< proto::terminal<proto::_>, proto::if_<IsArray<proto::_value>()> > {}; 

 typedef indexmaps::cuboid_domain<0> CD0;

 template<typename T> struct wrap_scalar { 
  typedef T value_type; 
  typedef CD0 domain_type;
  wrap_scalar( T const &) {}
  domain_type domain() const { return CD0();}
 };

 template<typename D1,typename D2> D1 _do (D1 const & a, D2 const & b) { 
   TRIQS_ARRAYS_DEBUG_CHECK((a==b), "Error : combining 2 arrays with different domains sizes"<<a<<" "<<b);
   return a;
  }
 template< typename D2> D2 _do (CD0 const & a, D2 const & b) { return b;} 
 template< typename D1> D1 _do (D1 const & a, CD0 const & b) { return a;} 

 template<typename ProtoTag, typename T1, typename T2> struct TypeAndDomain { 
  typedef typename T1::value_type V1; typedef typename T1::domain_type D1;
  typedef typename T2::value_type V2; typedef typename T2::domain_type D2;
  typedef BOOST_TYPEOF_TPL( V1() + V2()) value_type;
  static_assert ( (D1::rank==0) || (D2::rank==0) || (boost::is_same<D1,D2>::value), "Domain type mismatch");
  typedef typename mpl::if_c<(D1::rank!=0), D1,D2>::type domain_type;
  T1 const & a; T2 const & b;
  TypeAndDomain (T1 const & a_, T2 const & b_):a(a_),b(b_) {} 
  domain_type domain() const { return _do(a.domain(),b.domain()); }
 }; 

 struct ArrayGrammar : 
  proto::or_<
  proto::when< dC_ScalarGrammar,wrap_scalar<proto::_value>(proto::_value) >
  ,proto::when< BasicArrayTypeGrammar,proto::_value >
  ,proto::when< proto::plus <ArrayGrammar,ArrayGrammar>,TypeAndDomain<p_tag::plus,_left,_right > (_left,_right) >
  ,proto::when< proto::minus <ArrayGrammar,ArrayGrammar>,TypeAndDomain<p_tag::minus,_left,_right > (_left,_right)  >
  ,proto::when< proto::multiplies<ArrayGrammar,ArrayGrammar>,TypeAndDomain< p_tag::multiplies,_left,_right > (_left,_right)>
  ,proto::when< proto::divides<ArrayGrammar,ArrayGrammar>,TypeAndDomain< p_tag::divides,_left,_right > (_left,_right)>
  ,proto::when< proto::negate<ArrayGrammar >,ArrayGrammar(proto::_child0)  >
  > {};

 template<typename Expr> struct ArrayExpr;
 struct ArrayDomain : proto::domain<proto::generator<ArrayExpr>, ArrayGrammar> {};

 //   Evaluation context
 template<typename KeyType, typename ReturnType>
  struct ArrayEvalCtx : proto::callable_context< ArrayEvalCtx<KeyType,ReturnType> const > {
   typedef ReturnType result_type;
   KeyType const & key;
   ArrayEvalCtx(KeyType const & key_) : key(key_) {}
   template<typename T>// overrule just the terminals which have array interface.
    typename boost::enable_if< has_immutable_array_interface<T>, result_type >::type 
    operator ()(p_tag::terminal, T const & t) const { return t[key]; }
  };

 //   Expression
 template<typename Expr> struct ArrayExpr : Tag::expression, proto::extends<Expr, ArrayExpr<Expr>, ArrayDomain> { 
  typedef proto::extends<Expr, ArrayExpr<Expr>, ArrayDomain> base_type;
  typedef typename boost::result_of<ArrayGrammar(Expr) >::type _T;
  typedef typename _T::value_type value_type;
  typedef typename _T::domain_type domain_type;
  typedef typename domain_type::index_value_type key_type;
  typedef size_t index_type;
  static const int rank = _T::domain_type::rank;
  
  ArrayExpr( Expr const & expr = Expr() ) : base_type( expr ) {}
  domain_type domain() const { return ArrayGrammar()(*this).domain(); }
  value_type operator[] (key_type const &key) const { return proto::eval(*this, ArrayEvalCtx<key_type ,value_type> (key)); }
  value_type operator()(index_type const & i1) const { static_assert( (rank==1), "1 index expected"); return (*this)[mini_vector<size_t,1>(i1)]; }
  value_type operator()(index_type const & i1, index_type const & i2) const { static_assert((rank==2), "2 indices expected"); return (*this)[mini_vector<size_t,2>(i1,i2)]; }
  friend std::ostream &operator <<(std::ostream &sout, ArrayExpr<Expr> const &expr) { return proto::eval(expr, AlgebraPrintCtx (sout)); }
 };
}}

BOOST_PROTO_DEFINE_OPERATORS(expressions::array_algebra::IsArray, expressions::array_algebra::ArrayDomain);

template<typename Expr > array_view <typename Expr::value_type, Expr::domain_type::rank>
eval( Expr const & e) { return array<typename Expr::value_type, Expr::domain_type::rank>(e);}

}}//namespace triqs::arrays

#endif
