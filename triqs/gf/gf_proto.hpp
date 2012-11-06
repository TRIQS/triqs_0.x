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
#ifndef TRIQS_GF_PROTO_H
#define TRIQS_GF_PROTO_H
#include "./tools.hpp"

namespace triqs { namespace gf { 

 struct gf_no_mesh{}; // absence of mesh (for scalar)
 struct gf_no_shape{}; // absence of shape (for scalar)

 template< template< typename Expr> class scalar_identification_trait, template< typename Expr> class gf_identify, int gf_arity> 
  struct gf_proto_tools {  

   struct gf_scalar_grammar : proto::and_< proto::terminal<proto::_>, proto::if_<scalar_identification_trait<proto::_value>()> > {}; 

   struct gf_leaf_grammar: proto::and_< proto::terminal<proto::_>, proto::if_<gf_identify<proto::_value>()> > {}; 

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

   //------ computation of the data of the expression ----

   struct gf_get_data_tr {
    BOOST_PROTO_CALLABLE();
    template<typename Sig> struct result;
    template<typename This, typename F> struct result<This(F)> { typedef decltype( std::declval<const F>().data_view()) type; };
    template<typename F> typename result<gf_get_data_tr(F)>::type operator ()(F const &f) const { return f.data_view(); }
   };

   struct gf_data_tr : // evaluation of an expression  
    proto::or_<
    proto::when<gf_scalar_grammar, proto::_value>
    ,proto::when<gf_leaf_grammar, gf_get_data_tr(proto::_value) >
    ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::_default<gf_data_tr > >
    > {};

   //------ computation of the singularity of the expression ----

   struct gf_get_singularity_tr {
    BOOST_PROTO_CALLABLE();
    template<typename Sig> struct result;
    template<typename This, typename F> struct result<This(F)> { typedef decltype( std::declval<const F>().singularity_view()) type; };
    template<typename F> typename result<gf_get_singularity_tr(F)>::type operator ()(F const &f) const { return f.singularity_view(); }
   };

   struct gf_singularity_tr : // evaluation of an expression  
    proto::or_<
    proto::when<gf_scalar_grammar, proto::_value>
    ,proto::when<gf_leaf_grammar, gf_get_singularity_tr(proto::_value) >
    ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::_default<gf_singularity_tr > >
    > {};

   //------ computation of the mesh of the expression ----

   struct gf_combine_mesh_t { // a transform that combine the mesh of a node with the State
    BOOST_PROTO_CALLABLE();
    template<typename Sig> struct result;
    template<typename This, typename X, typename S> struct result<This(X,S)> {
     typedef typename tup::remove_const_and_ref<X>::type::mesh_t type;
     typedef typename tup::remove_const_and_ref<S>::type S2;
     static_assert((mpl::or_<boost::is_same<S2,gf_no_mesh>, boost::is_same<type,S2> >::value), "FATAL : two meshs of different type mixed in an expression");   
    };
    template<typename X> typename X::mesh_t operator ()(X const & x, gf_no_mesh ) const { return x.mesh();}
    template<typename X, typename M> typename X::mesh_t operator ()(X const & x, M const & m ) const { 
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
/*
   struct gf_combine_shape_t { 
    BOOST_PROTO_CALLABLE();
    template<typename Sig> struct result;
    template<typename This, typename X, typename S> struct result<This(X,S)> {
     typedef typename tup::remove_const_and_ref<X>::type::shape_t type;
     typedef typename tup::remove_const_and_ref<S>::type S2;
     static_assert((mpl::or_<boost::is_same<S2,gf_no_shape>, boost::is_same<type,S2> >::value), "FATAL : two shapes of different type mixed in an expression");   
    };
    template<typename X> typename X::shape_t operator ()(X const & x, gf_no_shape ) const { return x.shape();}
    template<typename X, typename M> typename X::shape_t operator ()(X const & x, M const & m ) const { 
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
*/
  };
}}

// I need to use a macro unfortunately here
// Indeed gf_expr must be in the same namespace as gf<DESC>, for proper ADL 
// and gf_expr can not be templated on anything else than Expr...
// And I want a separate expression domain for each kind of GF to detect errors early
// Usage is : TRIQS_GF_DEFINE_OPERATORS (DescriptorName, ScalarIdentificationTrait, GFIdentificationTrait)

#define TRIQS_GF_DEFINE_OPERATORS(DESC,ST,GFT) \
 template <typename Expr> struct gf_expr_##DESC;\
typedef tup::domain_with_copy<gf_proto_tools<ST,GFT,DESC::arity>::gf_grammar,gf_expr_##DESC> gf_expr_domain_##DESC;\
template<typename Expr> struct gf_expr_##DESC : DESC::tag, proto::extends<Expr, gf_expr_##DESC<Expr>, gf_expr_domain_##DESC> {\
 \
 typedef gf_proto_tools<ST,GFT,DESC::arity>::gf_mesh_t gf_mesh_t;\
 typedef gf_proto_tools<ST,GFT,DESC::arity>::gf_eval_t gf_eval_t;\
 typedef gf_proto_tools<ST,GFT,DESC::arity>::gf_data_tr gf_data_tr;\
 typedef gf_proto_tools<ST,GFT,DESC::arity>::gf_singularity_tr gf_singularity_tr;\
 \
 gf_expr_##DESC ( Expr const & expr = Expr() ) : proto::extends<Expr, gf_expr_##DESC <Expr>, gf_expr_domain_##DESC> (expr) {}\
 \
 typedef typename boost::result_of<gf_singularity_tr(Expr) >::type singularity_t;\
 singularity_t singularity_view()   const {return gf_singularity_tr() (*this);} \
 \
 typedef typename boost::result_of<gf_data_tr(Expr) >::type data_t;\
 data_t data_view()   const {return gf_data_tr() (*this);} \
 \
 typedef typename boost::result_of<gf_mesh_t(Expr,gf_no_mesh) >::type mesh_t;\
 mesh_t mesh()   const {return gf_mesh_t() (*this,gf_no_mesh());} \
 \
template<typename T> typename boost::result_of<gf_eval_t(Expr,bf::vector<T>)>::type \
 operator()(T const & x) const {return gf_eval_t()(*this,bf::make_vector(x));}\
 \
 template<typename T1,typename T2> typename boost::result_of<gf_eval_t(Expr,bf::vector<T1,T2>)>::type \
 operator()(T1 const & x1, T2 const & x2) const { return gf_eval_t()(*this,bf::make_vector(x1,x2));}\
 \
 friend std::ostream &operator <<(std::ostream &sout, gf_expr_##DESC <Expr> const &expr) { return tup::print_algebra(sout,expr);} \
};\
BOOST_PROTO_DEFINE_OPERATORS(GFT, gf_expr_domain_##DESC);
 

/*
 * typedef gf_proto_tools<ST,GFT,DESC::arity>::gf_shape_t gf_shape_t;\
 typedef typename boost::result_of<gf_shape_t(Expr,gf_no_shape) >::type shape_t;\
 shape_t shape() const {return gf_shape_t()(*this,gf_no_shape());} \
 \
*/ 
#endif