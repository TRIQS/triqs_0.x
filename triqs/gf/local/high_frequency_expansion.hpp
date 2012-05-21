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
#ifndef TRIQS_GF_LOCAL_TAIL_H
#define TRIQS_GF_LOCAL_TAIL_H
#include "matrix_valued_fnt_on_mesh.hpp"

namespace triqs { namespace gf { namespace local {

 class tail;         // the value class
 class tail_view;    // the view class
 template<typename E>  class tail_expr; // the proto expression of tail

 template<typename G> struct LocalTail : mpl::false_{};  // a boolean trait to identify the objects modelling the concept LocalTail
 template<> struct LocalTail<tail >     : mpl::true_{};
 template<> struct LocalTail<tail_view >: mpl::true_{};
 template<typename M> struct LocalTail<tail_expr<M> >: mpl::true_{};

 typedef std::complex<double> dcomplex;
 typedef tqa::matrix<std::complex<double> >  mv_dcomplex_type;

 // ---------------------- implementation --------------------------------

 /// A common implementation class
 template<bool IsView> class tail_impl : 
  public   mv_fnt_on_mesh<dcomplex, meshes::tail,IsView, void > {
   typedef mv_fnt_on_mesh<dcomplex, meshes::tail,IsView, void > B;
   public : 
   tail_impl():B() {} 
   tail_impl(size_t N1, size_t N2, meshes::tail const & m):B(N1,N2,m) {}
   tail_impl(tail_impl const & m, range R1, range R2) : B(m,R1,R2){} // slice constructor 
   tail_impl(tail const & g): B(g){}
   tail_impl(tail_view const & g): B(g){} 
   //tail_impl(tail_impl const & g): B(g){}
   //tail_impl(tail_impl<!IsView> const & g): B(g){} 

   typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
   typedef tail_view view_type;
   typedef tail      non_view_type;

   template<typename RHS> void operator=(RHS rhs) {B::operator=(rhs);}

   using B::operator();// import evaluation on a mesh point 

   typename B::mv_type operator() (int n)       { 
    if (n>this->mesh().order_max()) TRIQS_RUNTIME_ERROR<<" n > Max Order "; 
    if (n<this->mesh().order_min()) TRIQS_RUNTIME_ERROR<<" n < Min Order ";
    return this->data(range(), range(), n- this->mesh().order_min());
   }

   typename B::const_mv_type operator() (int n) const { 
    if (n>this->mesh().order_max()) TRIQS_RUNTIME_ERROR<<" n > Max Order "; 
    if (n<this->mesh().order_min())  return typename B::mv_type::non_view_type();
    return this->data(range(), range(), n- this->mesh().order_min());
   }

   operator domains::infty() const { return domains::infty();}
  };

 ///The View class of GF
 class tail_view : public tail_impl <true> { 
  typedef tail_impl <true>  B;
  tail_view(tail_view const & m, range R1, range R2) : B(m,R1,R2){} // slice constructor 
  friend class tail; //friend class tail_impl;
  public :
  tail_view(tail_view const & g): B(g){}
  tail_view(tail const & g): B(g){}

  using B::operator=; // import operator = from impl. class or the default = is synthetized and is the only one
  using B::operator(); // import all previously defined operator() for overloading
  TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,tail_view);

  /// Slice in orbital space
  const view_type slice(range R1, range R2) const { return view_type(*this,R1,R2); } 

  friend std::ostream & triqs_nvl_formal_print(std::ostream & out, tail_view const & x) { return out<<"tail_view";}
  friend std::ostream & operator << (std::ostream & out, tail_view const & x) { return out<<"tail_view";}
 };

 ///The regular class 
 class tail : public tail_impl <false> { 
  typedef tail_impl <false>  B;
  public : 
  tail():B() {} 
  tail(size_t N1, size_t N2, B::mesh_type const & m):B(N1,N2,m) {} // CHANGE THIS 
  tail(tail const & g): B(g){}
  tail(tail_view const & g): B(g){} 
  template<typename GfType> tail(GfType const & x): B() { *this = x;} // to maintain value semantics

  using B::operator=;
  using B::operator();
  TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,tail_view);

  /// Slice in orbital space
  const view_type slice(range R1, range R2) const { return view_type(*this,R1,R2); } 
 };

 // -------------------------------   Expression template for tail and view -----------------------

 // -----------  tail special multiplication --------------------

 template<typename L, typename R> struct tail_mul_lazy { // should have the same concept as tail_expr
  typedef typename tup::const_view_type_if_exists_else_type<L>::type L_type; L_type l; 
  typedef typename tup::const_view_type_if_exists_else_type<R>::type R_type; R_type r;
  tail_mul_lazy(L const & l_, R const & r_):l(l_),r(r_){}
  typedef meshes::tail mesh_type;
  meshes::tail domain () const { 
   int omin = l.mesh().order_min()+r.mesh().order_min(); 
   return meshes::tail(omin,omin + std::min(l.mesh().size(),r.mesh().size()));
  }
  mv_dcomplex_type operator ()(int n) const {
   // sum_{p}^n a_p b_{n-p}. p <= a.n_max, p >= a.n_min and n-p <=b.n_max and n-p >= b.n_min
   // hence p <= min ( a.n_max, n-b.n_min ) and p >= max ( a.n_min, n- b.n_max)  
   const int pmin = std::max(l.mesh().order_min(), n - r.mesh().order_max() );
   const int pmax = std::min(l.mesh().order_max(), n - r.mesh().order_min() );
   mv_dcomplex_type::non_view_type res;
   if (pmin<=pmax) {
    res=l(pmin)*r(n-pmin);// force resize
    for (int p = pmin + 1; p <= pmax; ++p)  { res += l(p) * r(n-p);}   
   }
   else { 
    BOOST_AUTO( l0, l(l.mesh().order_min())); BOOST_AUTO( r0, r(r.mesh().order_min()));
    res.resize(l0.shape()[0],r0.shape()[1]); 
    res()=0; 
   }
   return res;
  }
 };

 // -----------  

 template<typename T1, typename T2> 
  typename boost::enable_if< mpl::and_<LocalTail<T1>, LocalTail<T2> >, tail_mul_lazy<T1,T2> >::type
  operator* (T1 const & a, T2 const & b) { return tail_mul_lazy<T1,T2>(a,b); }

 // -----------  tail special inversion --------------------

 template<typename T> struct tail_inv_lazy { // should have the same concept as tail_expr
  typename tup::const_view_type_if_exists_else_type<T>::type t; 
  tail_inv_lazy(tail_view const & t_):t(t_){}

  struct internal_data { // implementing the pattern LazyPreCompute
   tail t_inv;
   template<typename Ta> internal_data(Ta const & t): 
    t_inv(t.shape()[0], t.shape()[1], meshes::tail(- t.mesh().order_min(),- t.mesh().order_min() + t.mesh().size())) {
     // compute the inverse
     // b_n = - a_0^{-1} * sum_{p=0}^{n-1} b_p a_{n-p} for n>0
     // b_0 = a_0^{-1}
     // b_min <= p <=b_max ;  a_min <= n-p <= a_max ---> n-a_max  <= p <= n-a_min 
     const int omin = t_inv.mesh().order_min(); const int omax = t_inv.mesh().order_max();
     t_inv(omin) = inverse(t(t.mesh().order_min()));
     for (int n=omin+1; n<=omax;n++) {
      const int pmin = std::max(omin, n - t.mesh().order_max() );
      const int pmax = std::min(n-1, n - t.mesh().order_min() );
      for (int p=pmin; p<=pmax; p++) t_inv(n) -= t_inv(p)*t(n-p); //n-p >= t_omin --> p <= n-t_omin 
      t_inv(n) *= t_inv(omin);
     }
    }
  };
  friend struct internal_data;
  mutable boost::shared_ptr<internal_data> _id;
  void activate() const { if (!_id) _id= boost::make_shared<internal_data>(t);}

  mv_dcomplex_type operator ()(int n) const { activate(); return _id->t_inv(n); }
 };

 // -----------  

 template<typename T> typename boost::enable_if< LocalTail<T>, tail_inv_lazy <T> >::type 
  inverse (T const & t) { return tail_inv_lazy<T>(t);}

 template<typename A, typename T> // anything / tail ---> anything * inverse(tail)
  typename boost::lazy_enable_if< LocalTail<T>, tup::type_of_mult<A, tail_inv_lazy <T> > >::type 
  operator/ (A const & a, T const & t) { return a * tail_inv_lazy<T>(t);}

 // -----------  expression  --------------------

 // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
 template <typename T> struct is_scalar_or_element : mpl::or_< tqa::is_matrix_expr<T>, tup::is_in_ZRC<T> > {};
 struct tail_scalar_grammar : proto::and_< proto::terminal<proto::_>, proto::if_<is_scalar_or_element<proto::_value>()> > {}; 

 struct tail_leaf_grammar: proto::and_< proto::terminal<proto::_>, proto::if_<LocalTail<proto::_value>()> > {}; 

 struct tail_grammar : 
  proto::or_<
  tail_scalar_grammar , tail_leaf_grammar
  , proto::plus      <tail_grammar ,tail_grammar>
  , proto::minus     <tail_grammar ,tail_grammar>
  , proto::multiplies<tail_grammar ,tail_scalar_grammar>
  , proto::multiplies<tail_scalar_grammar,tail_scalar_grammar>
  , proto::divides   <tail_grammar ,tail_scalar_grammar>
  , proto::negate    <tail_grammar >
  > {};

 struct tail_eval_scalar_t { // a transform that evaluates a scalar as a trivial expansion (only order 0)
  BOOST_PROTO_CALLABLE();
  template<typename Sig> struct result;
  template<typename This, typename S, typename AL> struct result<This(S,AL)> { typedef tail::value_type type;};
  template<typename S, typename AL> tail::value_type operator ()(S const & s, AL const & al) const 
  { int n = bf::at<mpl::int_<0> >(al); return (n==0 ? s : S()); }
 };

 struct tail_eval_t : // evaluation of an expression. same as for gf.  
  proto::or_<
  proto::when<tail_scalar_grammar, tail_eval_scalar_t(proto::_value, proto::_state)>
  ,proto::when<tail_leaf_grammar, tup::eval_fnt<1>(proto::_value, proto::_state) >
  ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::_default<tail_eval_t > >
  > {};

 //------ computation of the mesh of the expression. similar as for gf, except the run time check is more relax ----

 struct no_mesh{}; // absence of mesh (for scalar)

 struct tail_combine_mesh_t { // a transform that combine the mesh of a node with the State
  BOOST_PROTO_CALLABLE();
  template<typename Sig> struct result;
  template<typename This, typename X, typename S> struct result<This(X,S)> {
   typedef meshes::tail type;
   typedef typename tup::remove_const_and_ref<S>::type S2;
   static_assert((mpl::or_<boost::is_same<S2,no_mesh>, boost::is_same<meshes::tail,S2> >::value), "FATAL : two meshs of different type mixed in an expression");   
  };
  template<typename X> typename X::mesh_type operator ()(X const & x, no_mesh ) const { return x.mesh();}
  template<typename X, typename M> typename X::mesh_type operator ()(X const & x, M const & m ) const { 
   return meshes::tail( std::min (x.mesh().order_min(), m.order_min()), std::min (x.mesh().order_max(), m.order_max())); 
  }
 };

 struct tail_mesh_t : // the transform that computes the mesh recursively using a fold and a State initialized to no_mesh
  proto::or_<
  proto::when < tail_scalar_grammar, proto::_state >
  ,proto::when< tail_leaf_grammar, tail_combine_mesh_t (proto::_value, proto::_state) >
  ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::fold<proto::_, proto::_state, tail_mesh_t >() >
  > {};

 // grammar and domain and proto expression
 typedef tup::domain_with_copy<tail_grammar,tail_expr> tail_expr_domain;

 template<typename Expr> struct tail_expr : proto::extends<Expr, tail_expr<Expr>, tail_expr_domain> {
  tail_expr( Expr const & expr = Expr() ) : proto::extends<Expr, tail_expr<Expr>, tail_expr_domain> (expr) {}
  typedef typename boost::result_of<tail_mesh_t(Expr,no_mesh) >::type mesh_type;
  mesh_type mesh() const {return tail_mesh_t()(*this,no_mesh());} 
  template<typename T> typename boost::result_of<tail_eval_t(Expr,bf::vector<T>)>::type 
   operator()(T const & x) const {return tail_eval_t()(*this,bf::make_vector(x));}
  friend std::ostream &operator <<(std::ostream &sout, tail_expr<Expr> const &expr) { return tup::print_algebra(sout,expr);} 
 };

 BOOST_PROTO_DEFINE_OPERATORS(LocalTail, tail_expr_domain);

}}}
#endif
