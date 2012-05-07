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
#ifndef TRIQS_GF_LOCAL_GF_H
#define TRIQS_GF_LOCAL_GF_H
#include "matrix_valued_fnt_on_mesh.hpp"

namespace triqs { namespace gf { namespace local {
 
 namespace tag { struct is_any_gf{}; template<typename MeshType> struct is_gf: is_any_gf {}; struct is_tail{}; }
 template <typename G> struct LocalGf   : boost::is_base_of< tag::is_any_gf, G> {}; // identification trait
 template <typename G> struct is_tail : boost::is_base_of< tag::is_tail  , G> {}; // identification trait

 template<typename MeshType> class gf;
 template<typename MeshType> class gf_view;

 typedef gf     <meshes::tail> high_frequency_expansion;
 typedef gf_view<meshes::tail> high_frequency_expansion_view;

 // ---------------------- implementation --------------------------------

 // derive from this to add a cast operator to infty only for tail... 
 template<typename MeshType> struct _add_infty_cast : tag::is_gf<MeshType> {}; 
 template<> struct _add_infty_cast<meshes::tail> : tag::is_any_gf  { operator domains::infty() const { return domains::infty();} };
 //template<> struct _add_infty_cast<meshes::tail> : tag::is_tail  { operator domains::infty() const { return domains::infty();} };

 ///The regular class of GF
 template<typename MeshType> class gf : 
  public _add_infty_cast<MeshType>, 
  public   mv_fnt_on_mesh<MeshType,false, typename mpl::if_c<MeshType::domain_type::has_tail,high_frequency_expansion, void>::type > {
   typedef mv_fnt_on_mesh<MeshType,false, typename mpl::if_c<MeshType::domain_type::has_tail,high_frequency_expansion, void>::type > B;
   public : 
   gf():B() {} 
   gf (size_t N1, size_t N2, typename B::mesh_type const & m):B(N1,N2,m) {}
   gf(gf const & g): B(g){}
   gf(gf_view<MeshType> const & g): B(g){} // for clarity only, the next construct could handle every case...
   template<typename GfType> gf(GfType const & x): B(x,0) {}

   typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
   typedef gf_view<MeshType> view_type;
   typedef gf<MeshType>      non_view_type;

   using B::operator=; // or the default one will synthetised and will be the only one !
   using B::operator(); TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

   template<typename A, typename B> const view_type slice(A a, B b) const { return view_type(*this,a,b); } 
   template<typename A, typename B> friend view_type slice(gf const & g, A a, B b) { g.slice(g,a,b); } 
  };

 // ---------------------------------------------------------------------------------
 ///The View class of GF
 template<typename MeshType> class gf_view : 
  public _add_infty_cast<MeshType>,
  public   mv_fnt_on_mesh<MeshType,true, typename mpl::if_c<MeshType::domain_type::has_tail,high_frequency_expansion, void>::type > {
   typedef mv_fnt_on_mesh<MeshType,true, typename mpl::if_c<MeshType::domain_type::has_tail,high_frequency_expansion, void>::type > B;
   protected:
   template<typename A, typename C> gf_view (gf_view const & m, A a, C c) : B(m,a,c){} 
   friend class gf<MeshType>;
   public :
   gf_view(gf_view const & g): B(g){}
   gf_view(gf<MeshType> const & g): B(g){}

   typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
   typedef gf_view<MeshType> view_type;
   typedef gf<MeshType>      non_view_type;

   using B::operator=; // or the default one will synthetised and will be the only one !
   using B::operator(); TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

   template<typename A, typename B> const view_type slice(A a, B b) const { return view_type(*this,a,b); } 
   template<typename A, typename B> friend view_type slice(gf_view const & g, A a, B b) { return g.slice(a,b); } 
   
   friend std::ostream & triqs_nvl_formal_print(std::ostream & out, gf_view const & x) { return out<<"gf_view";}
   friend std::ostream & operator << (std::ostream & out, gf_view const & x) { return out<<"gf_view";}
  };

 // -------------------------------   Expression template for gf  --------------------------------------------------

 namespace gf_expr_detail { 
  // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
  template <typename T> struct is_scalar_or_element : mpl::or_< tqa::is_matrix_expr<T>, tup::is_in_ZRC<T> > {};

  struct ElementGrammar: proto::and_< proto::terminal<proto::_>, proto::if_<LocalGf<proto::_value>()> > {}; 
  struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<is_scalar_or_element<proto::_value>()> > {}; 

  struct gf_grammar : // the grammar of gf
   proto::or_<
   ScalarGrammar , ElementGrammar
   , proto::plus      <gf_grammar,gf_grammar>
   , proto::minus     <gf_grammar,gf_grammar>
   , proto::multiplies<gf_grammar,gf_grammar>
   , proto::divides   <gf_grammar,gf_grammar>
   , proto::negate    <gf_grammar >
   > {};

  struct eval_t : // evaluation of an expression  
   proto::or_<
   proto::when<ScalarGrammar, proto::_value>
   ,proto::when<ElementGrammar, tup::eval_fnt<1>(proto::_value, proto::_state) >
   ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::_default<eval_t > >
   > {};

  //------ computation of the mesh of the expression ----

  struct no_mesh{}; // absence of mesh (for scalar)

  struct combine_mesh_t { // a transform that combine the mesh of a node with the State
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename X, typename S> struct result<This(X,S)> {
    typedef typename tup::remove_const_and_ref<X>::type::mesh_type type;
    typedef typename tup::remove_const_and_ref<S>::type S2;
    static_assert((mpl::or_<boost::is_same<S2,no_mesh>, boost::is_same<type,S2> >::value), "FATAL : two meshs of different type mixed in an expression");   
   };
   template<typename X> typename X::mesh_type operator ()(X const & x, no_mesh ) const { return x.mesh();}
   template<typename X, typename M> typename X::mesh_type operator ()(X const & x, M const & m ) const { 
    //if (!(x.mesh() == m)) TRIQS_RUNTIME_ERROR << "Domain mismatch in expressions"; // run time check 
    return x.mesh();
   }
  };

  struct mesh_t : // the transform that computes the mesh recursively using a fold and a State initialized to no_mesh
   proto::or_<
   proto::when < ScalarGrammar, proto::_state >
   ,proto::when< ElementGrammar, combine_mesh_t (proto::_value, proto::_state) >
   ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::fold<proto::_, proto::_state, mesh_t >() >
   > {};

  // grammar and domain and proto expression
  template <typename Expr> struct gf_expr;
  typedef tup::domain_with_copy<gf_grammar,gf_expr> gf_expr_domain;

  template<typename Expr> struct gf_expr : proto::extends<Expr, gf_expr<Expr>, gf_expr_domain> {
   gf_expr( Expr const & expr = Expr() ) : proto::extends<Expr, gf_expr<Expr>, gf_expr_domain> (expr) {}
   typedef typename boost::result_of<mesh_t(Expr,no_mesh) >::type mesh_type;
   mesh_type mesh() const {return mesh_t()(*this,no_mesh());} 
   template<typename T> typename boost::result_of<eval_t(Expr,bf::vector<T>)>::type 
    operator()(T const & x) const {return eval_t()(*this,bf::make_vector(x));}
   friend std::ostream &operator <<(std::ostream &sout, gf_expr<Expr> const &expr) { return tup::print_algebra(sout,expr);} 
  };
 }
 template <typename Expr> struct LocalGf<gf_expr_detail::gf_expr< Expr > >:mpl::true_{};

 BOOST_PROTO_DEFINE_OPERATORS(LocalGf, gf_expr_detail::gf_expr_domain);
}}}
#endif
