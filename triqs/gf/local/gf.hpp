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
//#include "high_frequency_expansion.hpp"

namespace triqs { namespace gf { namespace local {
 
 namespace tag { struct is_any_gf{}; template<typename MeshType> struct is_gf: is_any_gf {}; struct is_tail{}; }
 template <typename G> struct is_gf   : boost::is_base_of< tag::is_any_gf, G> {}; // identification trait
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

   friend std::ostream & triqs_nvl_formal_print(std::ostream & out, gf_view const & x) { return out<<"gf_view";}
   friend std::ostream & operator << (std::ostream & out, gf_view const & x) { return out<<"gf_view";}

   typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
   typedef gf_view<MeshType> view_type;
   typedef gf<MeshType>      non_view_type;

   using B::operator=; // or the default one will synthetised and will be the only one !
   using B::operator(); TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

   template<typename A, typename B> const view_type slice(A a, B b) const { return view_type(*this,a,b); } 
   template<typename A, typename B> friend view_type slice(gf_view const & g, A a, B b) { return g.slice(a,b); } 
  };

 // -------------------------------   Expression template for gf  --------------------------------------------------

 namespace proto=boost::proto;

 // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
 template <typename T> struct is_scalar_or_element : mpl::or_< tqa::is_matrix_expr<T>, triqs::utility::proto::is_in_ZRC<T> > {};

 namespace gf_expr_temp { 

  struct ElementGrammar: proto::and_< proto::terminal<proto::_>, proto::if_<is_gf<proto::_value>()> > {}; 
  struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<is_scalar_or_element<proto::_value>()> > {}; 

  struct gf_grammar : 
   proto::or_<
   ScalarGrammar , ElementGrammar
   , proto::plus      <gf_grammar,gf_grammar>
   , proto::minus     <gf_grammar,gf_grammar>
   , proto::multiplies<gf_grammar,gf_grammar>
   , proto::divides   <gf_grammar,gf_grammar>
   , proto::negate    <gf_grammar >
   > {};

  struct eval_t : 
   proto::or_<
   proto::when<ScalarGrammar, proto::_value>
   ,proto::when<ElementGrammar, tup::eval_fnt<1>(proto::_value, proto::_state) >
   ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::_default<eval_t > >
   > {};

  struct no_mesh{ typedef void domain_type;}; 

//#define V1 
#ifdef V1
  struct get_mesh {
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename X> struct result<This(X)> { typedef typename boost::remove_reference<X>::type::mesh_type type;};
   template<typename X> typename X::mesh_type  operator ()(X const & x) const { return x.mesh();}
  };

  struct combine_mesh {
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename M1, typename M2> struct result<This(M1,M2)> {typedef M2 type;};
   template<typename This, typename M1> struct result<This(M1,no_mesh)> {typedef M1 type;};
   template<typename M> M operator ()(no_mesh const & m1, M const & m2) const { return m2;}
   template<typename M> M operator ()(M const & m1, no_mesh const & m2) const { return m1;}
   template<typename M1, typename M2> M1 operator ()(M1 const & m1, M2 const & m2) const { 
    static_assert((boost::is_same<M1,M2>::value), "FATAL : two meshes of different type mixed in an expression");
    return m1;
   }
  };

  struct dom_t : 
   proto::or_<
   proto::when< ScalarGrammar, no_mesh() >
   ,proto::when< ElementGrammar, get_mesh(proto::_value) >
   ,proto::when< proto::plus <dom_t,dom_t>,       combine_mesh (dom_t(proto::_left), dom_t( proto::_right)) >
   ,proto::when< proto::minus <dom_t,dom_t>,      combine_mesh (dom_t(proto::_left), dom_t( proto::_right)) >
   ,proto::when< proto::multiplies <dom_t,dom_t>, combine_mesh (dom_t(proto::_left), dom_t( proto::_right)) >
   ,proto::when< proto::divides <dom_t,dom_t>,    combine_mesh (dom_t(proto::_left), dom_t( proto::_right)) >
   ,proto::when< proto::unary_expr<proto::_,dom_t >,  dom_t(proto::_left) >
   > {};

#else

  struct combine_mesh {
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename X, typename S> struct result<This(X,S)> {typedef typename tup::remove_const_and_ref<X>::type::mesh_type type;};
   template<typename X> typename X::mesh_type operator ()(X const & x, no_mesh ) const { return x.mesh();}
   template<typename X, typename M> typename X::mesh_type operator ()(X const & x, M const & m ) const { 
    //static_assert((boost::is_same<M1,M2>::value), "FATAL : two meshes of different type mixed in an expression");
    // run time check 
    return x.mesh();
   }
  };
  struct dom_t : 
   proto::or_<
   proto::when < ScalarGrammar, proto::_state >
   ,proto::when< ElementGrammar, combine_mesh (proto::_value, proto::_state) >
   ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::fold<proto::_, proto::_state, dom_t >() >
   > {};

#endif

  // grammar and domain and proto expression
  template <typename Expr> struct gf_expr;
  typedef tup::domain<gf_grammar,gf_expr,true> gf_expr_domain;

  typedef eval_t eval1;
  template<typename Expr> struct gf_expr : boost::proto::extends<Expr, gf_expr<Expr>, gf_expr_domain>{
   gf_expr( Expr const & expr = Expr() ) : boost::proto::extends<Expr, gf_expr<Expr>, gf_expr_domain>  ( expr ) {}
   typedef typename boost::result_of<dom_t(Expr,no_mesh) >::type mesh_type;
   typedef typename mesh_type::domain_type domain_type;
   mesh_type mesh() const { return dom_t()(*this,no_mesh()); } 
   template<typename T> 
    typename boost::result_of<eval1(Expr,bf::vector<T>) >::type operator() (T const & x) const {
     static_assert(mesh_type::domain_type::has_tail || !(boost::is_same<T,meshes::tail>::value), "Error");
     return eval1()(*this, bf::make_vector(x));
    }
   typedef tqa::array_view<std::string,2> indices_type;
   indices_type indices() const { return indices_type::non_view_type();}
   friend std::ostream &operator <<(std::ostream &sout, gf_expr<Expr> const &expr) { return boost::proto::eval(expr, tup::algebra_print_ctx (sout)); }
  };
 }
 template <typename Expr> struct is_gf<gf_expr_temp::gf_expr< Expr > >:mpl::true_{};

 BOOST_PROTO_DEFINE_OPERATORS(is_gf, gf_expr_temp::gf_expr_domain);
}}}
#endif
