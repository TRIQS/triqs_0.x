/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_GF_GFPROTO_H
#define TRIQS_GF_GFPROTO_H

/*This header has intentionnally no namespace, it is intended to be included *in* a namespace gf::local, etc....
 * One must define 
 * scalar_identification_trait: to identify the scalar
 *  gf_identification_trait: to identify the gf
 *  GfArity 
 */
struct gf_scalar_grammar : proto::and_< proto::terminal<proto::_>, proto::if_<scalar_identification_trait<proto::_value>()> > {}; 

struct gf_leaf_grammar: proto::and_< proto::terminal<proto::_>, proto::if_<gf_identification_trait<proto::_value>()> > {}; 

struct gf_grammar : // the grammar of gf
 proto::or_< gf_scalar_grammar , gf_leaf_grammar
 , proto::plus      <gf_grammar,gf_grammar>
 , proto::minus     <gf_grammar,gf_grammar>
 , proto::multiplies<gf_grammar,gf_grammar>
 , proto::divides   <gf_grammar,gf_grammar>
 , proto::negate    <gf_grammar >
 > {};

struct gf_eval_t : // evaluation of an expression  
 proto::or_<
 proto::when<gf_scalar_grammar, proto::_value>
 ,proto::when<gf_leaf_grammar, tup::eval_fnt<gf_arity>(proto::_value, proto::_state) >
 ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::_default<gf_eval_t > >
 > {};

//------ computation of the mesh of the expression ----

struct gf_no_mesh{}; // absence of mesh (for scalar)

struct gf_combine_mesh_t { // a transform that combine the mesh of a node with the State
 BOOST_PROTO_CALLABLE();
 template<typename Sig> struct result;
 template<typename This, typename X, typename S> struct result<This(X,S)> {
  typedef typename tup::remove_const_and_ref<X>::type::mesh_type type;
  typedef typename tup::remove_const_and_ref<S>::type S2;
  static_assert((mpl::or_<boost::is_same<S2,gf_no_mesh>, boost::is_same<type,S2> >::value), "FATAL : two meshs of different type mixed in an expression");   
 };
 template<typename X> typename X::mesh_type operator ()(X const & x, gf_no_mesh ) const { return x.mesh();}
 template<typename X, typename M> typename X::mesh_type operator ()(X const & x, M const & m ) const { 
  if (!(x.mesh() == m)) TRIQS_RUNTIME_ERROR << "Domain mismatch in expressions"; // run time check 
  return x.mesh();
 }
};

struct gf_mesh_t : // the transform that computes the mesh recursively using a fold and a State initialized to gf_no_mesh
 proto::or_<
 proto::when < gf_scalar_grammar, proto::_state >
 ,proto::when< gf_leaf_grammar, gf_combine_mesh_t (proto::_value, proto::_state) >
 ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::fold<proto::_, proto::_state, gf_mesh_t >() >
 > {};

//------ computation of the shape of the expression ----

struct gf_no_shape{}; // absence of shape (for scalar)

struct gf_combine_shape_t { 
 BOOST_PROTO_CALLABLE();
 template<typename Sig> struct result;
 template<typename This, typename X, typename S> struct result<This(X,S)> {
  typedef typename tup::remove_const_and_ref<X>::type::shape_type type;
  typedef typename tup::remove_const_and_ref<S>::type S2;
  static_assert((mpl::or_<boost::is_same<S2,gf_no_shape>, boost::is_same<type,S2> >::value), "FATAL : two shapes of different type mixed in an expression");   
 };
 template<typename X> typename X::shape_type operator ()(X const & x, gf_no_shape ) const { return x.shape();}
 template<typename X, typename M> typename X::shape_type operator ()(X const & x, M const & m ) const { 
  if (!(x.shape() == m)) TRIQS_RUNTIME_ERROR << "Domain mismatch in expressions"; // run time check 
  return x.shape();
 }
};

struct gf_shape_t : // the transform that computes the shape recursively using a fold and a State initialized to gf_no_shape
 proto::or_<
 proto::when < gf_scalar_grammar, proto::_state >
 ,proto::when< gf_leaf_grammar, gf_combine_shape_t (proto::_value, proto::_state) >
 ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::fold<proto::_, proto::_state, gf_shape_t >() >
 > {};

//------  grammar and domain and proto expression ----------------------

template <typename Expr> struct gf_expr;
typedef tup::domain_with_copy<gf_grammar,gf_expr> gf_expr_domain;

template<typename Expr> struct gf_expr : proto::extends<Expr, gf_expr<Expr>, gf_expr_domain> {

 gf_expr( Expr const & expr = Expr() ) : proto::extends<Expr, gf_expr<Expr>, gf_expr_domain> (expr) {}

 typedef typename boost::result_of<gf_mesh_t(Expr,gf_no_mesh) >::type mesh_type;
 mesh_type mesh()   const {return gf_mesh_t() (*this,gf_no_mesh());} 

 typedef typename boost::result_of<gf_shape_t(Expr,gf_no_shape) >::type shape_type;
 shape_type shape() const {return gf_shape_t()(*this,gf_no_shape());} 

 template<typename T> typename boost::result_of<gf_eval_t(Expr,bf::vector<T>)>::type 
  operator()(T const & x) const {return gf_eval_t()(*this,bf::make_vector(x));}

 friend std::ostream &operator <<(std::ostream &sout, gf_expr<Expr> const &expr) { return tup::print_algebra(sout,expr);} 
};

BOOST_PROTO_DEFINE_OPERATORS(LocalGf, gf_expr_domain);

#endif
