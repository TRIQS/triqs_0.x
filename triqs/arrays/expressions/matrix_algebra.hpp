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
#ifndef TRIQS_ARRAYS_EXPRESSION_MATRIX_ALGEBRA_H
#define TRIQS_ARRAYS_EXPRESSION_MATRIX_ALGEBRA_H
#include <triqs/utility/proto/tools.hpp>
#include "../matrix.hpp"
#include "../linalg/matmul.hpp"
#include "../linalg/mat_vec_mul.hpp"
#include "../linalg/inverse.hpp"

namespace triqs { namespace arrays { 

 namespace bf = boost::fusion; namespace tup = triqs::utility::proto; namespace mpl = boost::mpl; namespace proto = boost::proto;

 template<typename M1, typename M2> // matrix * matrix
  typename boost::enable_if< mpl::and_<is_matrix_expr<M1>, is_matrix_expr<M2> >, matmul_lazy<M1,M2> >::type
  operator* (M1 const & a, M2 const & b) { return matmul_lazy<M1,M2>(a,b); }

 template<typename M, typename V> // matrix * vector
  typename boost::enable_if< mpl::and_<is_matrix_expr<M>, is_vector_expr<V> >, mat_vec_mul_lazy<M,V> >::type
  operator* (M const & m, V const & v) { return mat_vec_mul_lazy<M,V>(m,v); }

 template<typename A, typename M> // anything / matrix ---> anything * inverse(matrix)
  typename boost::lazy_enable_if< is_matrix_expr<M>, tup::type_of_mult<A, inverse_lazy <M> > >::type
  operator/ (A const & a, M const & m) { return a * inverse(m);}

 template<typename Expr> struct matrix_expr;
 template<typename Expr> struct vector_expr;

 namespace detail_expr_matvec{ 

  using boost::enable_if; using boost::enable_if_c; using proto::_left; using proto::_right; 
  using proto::_value; namespace p_tag= proto::tag; using boost::remove_reference;

  typedef indexmaps::cuboid_domain<2> matrix_domain_type;

  struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<tup::is_in_ZRC<proto::_value>()> > {}; 
  struct BasicVectorTypeGrammar :  proto::and_< proto::terminal<proto::_>, proto::if_<is_vector_expr<proto::_value>()> > {}; 
  struct BasicMatrixTypeGrammar :  proto::and_< proto::terminal<proto::_>, proto::if_<is_matrix_expr<proto::_value>()> > {}; 

  struct MatrixGrammar : 
   proto::or_<
   ScalarGrammar // scalar in injected in the algebra 
   , BasicMatrixTypeGrammar
   , proto::plus      <MatrixGrammar,MatrixGrammar>
   , proto::minus     <MatrixGrammar,MatrixGrammar>
   , proto::multiplies<ScalarGrammar,MatrixGrammar>
   , proto::multiplies<MatrixGrammar,ScalarGrammar>
   , proto::divides   <MatrixGrammar,ScalarGrammar>
   , proto::negate    <MatrixGrammar >
   > {};

  struct VectorGrammar : 
   proto::or_<
   BasicVectorTypeGrammar
   , proto::plus      <VectorGrammar,VectorGrammar>
   , proto::minus     <VectorGrammar,VectorGrammar>
   , proto::multiplies<ScalarGrammar,VectorGrammar>
   , proto::multiplies<VectorGrammar,ScalarGrammar>
   , proto::divides   <VectorGrammar,ScalarGrammar>
   , proto::negate    <VectorGrammar >
   > {};

  struct eval_scalar { // a transform that evaluates a scalar as an IDENTITY MATRIX !
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename S, typename AL> struct result<This(S,AL)> { typedef typename remove_reference<S>::type type;};
   template<typename S, typename AL> 
    typename result<eval_scalar(S,AL)>::type operator ()(S const & s, AL const & al) const 
    { BOOST_AUTO( key, bf::at<mpl::int_<0> >(al)); return (key[0]==key[1] ? s : S()); }
  };

  struct eval_mv { // a transform that evaluates a matrix or a vector on a key contains at head of a bf::vector
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename M, typename AL> struct result<This(M,AL)> { typedef typename remove_reference<M>::type::value_type type;};
   template<typename M, typename AL> 
    typename result<eval_mv(M,AL)>::type operator ()(M const & m, AL const & al) const { return m[bf::at<mpl::int_<0> >(al)]; }
  };

  struct eval_t;
  struct eval_t_cases { template <typename TAG> struct case_: proto::not_<proto::_>{}; };
  template<> struct eval_t_cases::case_<proto::tag::terminal> :  
   proto::or_<
   proto::when<ScalarGrammar, eval_scalar(proto::_value,proto::_state)>
   ,proto::when<BasicMatrixTypeGrammar, eval_mv (proto::_value, proto::_state) >
   ,proto::when<BasicVectorTypeGrammar, eval_mv (proto::_value, proto::_state) >
   > {};
  template<> struct eval_t_cases::case_<proto::tag::plus>  : proto::when< proto::plus <eval_t,eval_t>,  proto::_default<eval_t> > {};
  template<> struct eval_t_cases::case_<proto::tag::minus> : proto::when< proto::minus <eval_t,eval_t>, proto::_default<eval_t> > {};
  template<> struct eval_t_cases::case_<proto::tag::multiplies> :  
   proto::or_<
   proto::when< proto::multiplies<ScalarGrammar,eval_t>,  tup::multiplies_t (_value(_left),eval_t(_right)) > 
   //proto::when< proto::multiplies<ScalarGrammar,eval_t>,  proto::_default<proto::multiplies<_value(_left),eval_t(_right)> >  > 
   ,proto::when< proto::multiplies<eval_t,ScalarGrammar>, tup::multiplies_t (_value(_right),eval_t(_left))  > 
   > {};
  template<> struct eval_t_cases::case_<proto::tag::divides>:proto::when< proto::divides<eval_t,ScalarGrammar>, tup::divides_t(eval_t(_left),_value(_right))>{};
  template<> struct eval_t_cases::case_<proto::tag::negate> : proto::when< proto::negate<eval_t>, eval_t(_left)> {};
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
   ,proto::when< BasicMatrixTypeGrammar, get_domain(proto::_value) >
   ,proto::when< BasicVectorTypeGrammar, get_domain(proto::_value) >
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
  struct MatrixDomain : proto::domain<proto::generator<matrix_expr>, MatrixGrammar> {
   template< typename T > struct as_child : proto_base_domain::as_expr< T > {};
  };

  struct VectorDomain : proto::domain<proto::generator<vector_expr>, VectorGrammar> {
   template< typename T > struct as_child : proto_base_domain::as_expr< T > {};
  };

  //For matrices and vectors, special treatment : the regular classes are replaced by the corresponding const view
  template< typename T, typename Opt> struct MatrixDomain::as_child< matrix<T,Opt> > : 
   MatrixDomain::proto_base_domain::template as_expr< const matrix_view<T,Opt> >{};

  template< typename T, typename Opt> struct MatrixDomain::as_child< const matrix<T,Opt> > : 
   MatrixDomain::proto_base_domain::template as_expr< const matrix_view<T,Opt> >{};

  template< typename T, typename Opt> struct VectorDomain::as_child< vector<T,Opt> > : 
   VectorDomain::proto_base_domain::template as_expr< const vector_view<T,Opt> >{};

  template< typename T, typename Opt> struct VectorDomain::as_child< const vector<T,Opt> > : 
   VectorDomain::proto_base_domain::template as_expr< const vector_view<T,Opt> >{};
 }

 //   Expression for matrices
 template<typename Expr> struct matrix_expr : proto::extends<Expr, matrix_expr<Expr>, detail_expr_matvec::MatrixDomain>, Tag::expression { 
  matrix_expr( Expr const & expr = Expr() ) : proto::extends<Expr, matrix_expr<Expr>, detail_expr_matvec::MatrixDomain> ( expr ) {} 

  typedef size_t index_type;
  typedef typename boost::remove_reference<typename boost::result_of<detail_expr_matvec::dom_t(matrix_expr) >::type>::type domain_type;
  static const int rank = domain_type::rank; static_assert((rank==2)||(rank==0), "2 indices expected"); 
  typedef typename domain_type::index_value_type key_type; 
  typedef typename boost::remove_reference<typename boost::result_of<detail_expr_matvec::eval_t(matrix_expr,bf::vector<key_type>) >::type>::type value_type;

  domain_type domain() const { return detail_expr_matvec::dom_t()(*this); }
  mini_vector<size_t,2> shape() const { return this->domain().lengths();} 
  size_t dim0() const { return this->domain().lengths()[0];} 
  size_t dim1() const { return this->domain().lengths()[1];} 

  value_type operator[] (key_type const &key) const { return detail_expr_matvec::eval_t()(*this, bf::make_vector(key)); }
  value_type operator()(index_type const & i1, index_type const & i2) const { return (*this)[mini_vector<size_t,2>(i1,i2)]; }

  friend std::ostream &operator <<(std::ostream &sout, matrix_expr<Expr> const &expr){return proto::eval(expr,tup::algebra_print_ctx (sout));}
 };

 //   Expression for vectors
 template<typename Expr> struct vector_expr : proto::extends<Expr, vector_expr<Expr>, detail_expr_matvec::VectorDomain>, Tag::expression { 
  vector_expr( Expr const & expr = Expr() ) : proto::extends<Expr, vector_expr<Expr>, detail_expr_matvec::VectorDomain> ( expr ) {} 

  typedef size_t index_type;
  typedef typename boost::remove_reference<typename boost::result_of<detail_expr_matvec::dom_t(vector_expr) >::type>::type domain_type;
  static const int rank = domain_type::rank; static_assert((rank==1)||(rank==0), "1 indice expected"); 
  typedef typename domain_type::index_value_type key_type; 
  typedef typename boost::remove_reference<typename boost::result_of<detail_expr_matvec::eval_t(vector_expr,bf::vector<key_type>) >::type>::type value_type;

  domain_type domain() const { return detail_expr_matvec::dom_t()(*this); }
  mini_vector<size_t,1> shape() const { return this->domain().lengths();} 
  size_t dim0() const { return this->domain().lengths()[0];} 

  value_type operator[] (key_type const &key) const { return detail_expr_matvec::eval_t()(*this, bf::make_vector(key)); }
  value_type operator()(index_type const & i1) const { return (*this)[mini_vector<size_t,1>(i1)]; }

  friend std::ostream &operator <<(std::ostream &sout, vector_expr<Expr> const &expr){return proto::eval(expr,tup::algebra_print_ctx (sout));}
 };

 template<typename Expr> struct is_matrix_expr<matrix_expr<Expr> > : 
  mpl::not_<boost::proto::matches<Expr,detail_expr_matvec::ScalarGrammar> >{}; // every expression but a single scalar
 template<typename Expr> struct is_vector_expr<vector_expr<Expr> > : mpl::true_{};

 BOOST_PROTO_DEFINE_OPERATORS(is_matrix_expr, detail_expr_matvec::MatrixDomain);
 BOOST_PROTO_DEFINE_OPERATORS(is_vector_expr, detail_expr_matvec::VectorDomain);

 template<typename Expr > matrix_view <typename Expr::value_type> make_matrix( Expr const & e) { return matrix<typename Expr::value_type>(e);}
 template<typename Expr > vector_view <typename Expr::value_type> make_vector( Expr const & e) { return vector<typename Expr::value_type>(e);}

}}//namespace triqs::arrays 

#endif

