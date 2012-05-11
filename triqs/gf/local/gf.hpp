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
#include "high_frequency_expansion.hpp"

namespace triqs { namespace gf { namespace local {
 
 template<typename MeshType> class gf;         // the value class
 template<typename MeshType> class gf_view;    // the view class
 template<typename E>        class gf_expr; // the proto expression of gf

 template<typename G> struct LocalGf : mpl::false_{};  // a boolean trait to identify the objects modelling the concept LocalGf
 template<typename M> struct LocalGf<gf<M> >     : mpl::true_{};
 template<typename M> struct LocalGf<gf_view<M> >: mpl::true_{};
 template<typename M> struct LocalGf<gf_expr<M> >: mpl::true_{};

 // ---------------------- implementation --------------------------------

 // a domain_type ---> what is the type of the tail. All function have a tail, except Legendre
 template<typename M, bool V> struct gf_impl1 { typedef mv_fnt_on_mesh<dcomplex, M,V, tail> type;};
 //template<> struct tail_type_from_domain<domains::legendre> { typedef void type;};
 //specialize for matsubara time for double

 ///The regular class of GF
 template<typename MeshType> class gf : 
  public   gf_impl1<MeshType,false>::type  {
   typedef typename gf_impl1<MeshType,false>::type  B;
   public : 
   gf():B() {} 
   gf(size_t N1, size_t N2, MeshType const & m):B(N1,N2,m) {}
   gf(gf const & g): B(g){}
   gf(gf_view<MeshType> const & g): B(g){} 
   template<typename GfType> gf(GfType const & x): B() { *this = x;} // to maintain value semantics

   typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
   typedef gf_view<MeshType> view_type;
   typedef gf<MeshType>      non_view_type;

   template<typename RHS> void operator=(RHS rhs) {B::operator=(rhs);}

   using B::operator();// import evaluation on a mesh point 

   typedef typename MeshType::domain_type::point_type arg0_type;
   typename B::mv_type       operator() (arg0_type const & i)       { return this->mesh().interpolate(*this,i);}
   typename B::const_mv_type operator() (arg0_type const & i) const { return this->mesh().interpolate(*this,i);}

   TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

   /// Slice in orbital space
   const view_type slice(range R1, range R2) const { return view_type(*this,R1,R2); } 

   /// Slice of the mesh
   view_type       slice_mesh(typename MeshType::slice_arg_type arg)       { return view_type (*this, arg);}
   view_type const slice_mesh(typename MeshType::slice_arg_type arg) const { return view_type (*this, arg);}
  };

 // ---------------------------------------------------------------------------------
 ///The View class of GF
 template<typename MeshType> class gf_view : 
  public   gf_impl1<MeshType,true>::type  {
   typedef typename gf_impl1<MeshType,true>::type B;
   protected:
   gf_view (gf_view const & m, range R1, range R2) : B(m,R1,R2){} // slice constructor 
   friend class gf<MeshType>;
   public :
   gf_view(gf_view const & g): B(g){}
   gf_view(gf<MeshType> const & g): B(g){}

   typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
   typedef gf_view<MeshType> view_type;
   typedef gf<MeshType>      non_view_type;

   template<typename RHS> void operator=(RHS rhs) {B::operator=(rhs);}

   using B::operator();// import evaluation on a mesh point 

   typedef typename MeshType::domain_type::point_type arg0_type;
   typename B::mv_type       operator() (arg0_type const & i)       { return this->mesh().interpolate(*this,i);} // MUST BE different ! can not affect to tail !
   typename B::const_mv_type operator() (arg0_type const & i) const { return this->mesh().interpolate(*this,i);}

   TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

   /// Slice in orbital space
   const view_type slice(range R1, range R2) const { return view_type(*this,R1,R2); } 

   /// Slice of the mesh
   view_type       slice_mesh(typename MeshType::slice_arg_type arg)       { return view_type (*this, arg);}
   view_type const slice_mesh(typename MeshType::slice_arg_type arg) const { return view_type (*this, arg);}

   friend std::ostream & triqs_nvl_formal_print(std::ostream & out, gf_view const & x) { return out<<"gf_view";}
   friend std::ostream & operator << (std::ostream & out, gf_view const & x) { return out<<"gf_view";}
  };

 // -------------------------------   Expression template for gf  --------------------------------------------------

 // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
 //template <typename T> struct is_scalar_or_element : mpl::or_< tqa::is_matrix_expr<T>, tup::is_in_ZRC<T> > {};
 struct gf_scalar_grammar : proto::and_< proto::terminal<proto::_>, proto::if_<is_scalar_or_element<proto::_value>()> > {}; 

 struct gf_leaf_grammar: proto::and_< proto::terminal<proto::_>, proto::if_<LocalGf<proto::_value>()> > {}; 

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
  ,proto::when<gf_leaf_grammar, tup::eval_fnt<1>(proto::_value, proto::_state) >
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
   static_assert((mpl::or_<boost::is_same<S2,no_mesh>, boost::is_same<type,S2> >::value), "FATAL : two meshs of different type mixed in an expression");   
  };
  template<typename X> typename X::mesh_type operator ()(X const & x, no_mesh ) const { return x.mesh();}
  template<typename X, typename M> typename X::mesh_type operator ()(X const & x, M const & m ) const { 
   if (!(x.mesh() == m)) TRIQS_RUNTIME_ERROR << "Domain mismatch in expressions"; // run time check 
   return x.mesh();
  }
 };

 struct gf_mesh_t : // the transform that computes the mesh recursively using a fold and a State initialized to no_mesh
  proto::or_<
  proto::when < gf_scalar_grammar, proto::_state >
  ,proto::when< gf_leaf_grammar, gf_combine_mesh_t (proto::_value, proto::_state) >
  ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::fold<proto::_, proto::_state, gf_mesh_t >() >
  > {};

 // grammar and domain and proto expression
 template <typename Expr> struct gf_expr;
 typedef tup::domain_with_copy<gf_grammar,gf_expr> gf_expr_domain;

 template<typename Expr> struct gf_expr : proto::extends<Expr, gf_expr<Expr>, gf_expr_domain> {
  gf_expr( Expr const & expr = Expr() ) : proto::extends<Expr, gf_expr<Expr>, gf_expr_domain> (expr) {}
  typedef typename boost::result_of<gf_mesh_t(Expr,no_mesh) >::type mesh_type;
  mesh_type mesh() const {return gf_mesh_t()(*this,no_mesh());} 
  template<typename T> typename boost::result_of<gf_eval_t(Expr,bf::vector<T>)>::type 
   operator()(T const & x) const {return gf_eval_t()(*this,bf::make_vector(x));}
  friend std::ostream &operator <<(std::ostream &sout, gf_expr<Expr> const &expr) { return tup::print_algebra(sout,expr);} 
 };

 BOOST_PROTO_DEFINE_OPERATORS(LocalGf, gf_expr_domain);
}}}
#endif
