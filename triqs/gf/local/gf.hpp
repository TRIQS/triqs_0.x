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
#include "tail.hpp"

namespace triqs { namespace gf { namespace local {
 
 template<typename MeshType> class gf;         // the value class
 template<typename MeshType> class gf_view;    // the view class
 template<typename E>        class gf_expr; // the proto expression of gf

 template<typename G> struct LocalGf : mpl::false_{};  // a boolean trait to identify the objects modelling the concept LocalGf
 template<typename M> struct LocalGf<gf<M> >     : mpl::true_{};
 template<typename M> struct LocalGf<gf_view<M> >: mpl::true_{};
 template<typename M> struct LocalGf<gf_expr<M> >: mpl::true_{};

 // ---------------------- implementation --------------------------------

 // traits to compute ...
 // special case : gf in time
 template<typename M> struct value_type_from_mesh_type { typedef dcomplex type;};
 template<typename M> struct tail_type_from_mesh_type { typedef tail type;};

 ///The regular class of GF
 template<typename MeshType,bool IsView> class gf_impl {
   public : 
   typedef MeshType mesh_type;
   typedef typename mesh_type::domain_type domain_type;
   static const bool has_tail = mesh_type::domain_type::has_tail;

   typedef typename value_type_from_mesh_type<mesh_type>::type value_type;
   typedef arrays::Option::Fortran storage_order;
   typedef arrays::array      <value_type,3,storage_order>                       data_non_view_type;
   typedef arrays::array_view <value_type,3,storage_order>                       data_view_type;
   typedef typename mpl::if_c<IsView, data_view_type, data_non_view_type>::type  data_type;

   typedef typename tail_type_from_mesh_type<mesh_type>::type                     tail_non_view_type; 
   typedef typename tail_non_view_type::view_type                                 tail_view_type; 
   typedef typename mpl::if_c<IsView, tail_view_type, tail_non_view_type>::type   tail_type;
  
   mesh_type const & mesh() const { return mymesh;}
   domain_type const & domain() const { return mymesh.domain();}
   //data_view_type data_view()             { return data;} 
   //const data_view_type data_view() const { return data;}

   typedef tqa::mini_vector<size_t,2> shape_type;
   shape_type shape() const { return shape_type(data.shape()[0], data.shape()[1]);} 

   protected:
   mesh_type mymesh;
   data_type data;
   tail_type tail;
 
   gf_impl() {} // all arrays of zero size (empty)
   gf_impl(size_t N1, size_t N2, MeshType const & m, int tail_order_max): mymesh(m), data(N1,N2,m.size()), tail(N1,N2,-1,tail_order_max){ data()=0;}
   gf_impl (mesh_type const & m, data_view_type const & dat,tail_view_type const & t) : mymesh(m), data(dat), tail(t){}
   
   gf_impl(gf_impl const & x)                  : mymesh(x.mymesh), data(x.data), tail(x.tail){}// for clarity, would be synthetized anyway
   gf_impl(gf_impl<MeshType,!IsView> const & x): mymesh(x.mymesh), data(x.data), tail(x.tail){} 
   
   // orbital slice constructor. Only for view.
   template<bool V> gf_impl (gf_impl<MeshType,V> const & m, range r1, range r2) : 
     mymesh(m.mesh()),data(m.data(r1,r2, range())),tail(m.tail.slice(r1,r2)){ static_assert(IsView, "Internal Error");}
   
   // mesh slice constructor. Only for view.
   template<bool V> gf_impl (gf_impl<MeshType,V> const & m, typename MeshType::slice_arg_type arg) : 
     mymesh(m.mesh().slice(arg)),data(m.data),tail(m.tail){static_assert(IsView, "Internal Error");}

   friend class gf_impl<MeshType,!IsView>;

   // access to the data . Beware, we view it as a *matrix* NOT an array...
   tqa::matrix_view <value_type, storage_order> operator[](size_t u) { return data(range(),range(),u);}

   public:
   
   void operator = (gf_impl const & rhs) { mymesh = rhs.mymesh; data = rhs.data; tail = rhs.tail; }// for clarity, would be synthetized anyway
   void operator = (gf_impl<MeshType,!IsView> const & rhs) { mymesh = rhs.mymesh; data = rhs.data; tail = rhs.tail; }
   template<typename RHS> void operator = (RHS const & rhs) { // the general version  
    mymesh = rhs.mesh(); 
    BOOST_AUTO(r0 , rhs(mesh()[0])); // first point computed first to compute the shape of the result... 
    tqa::resize_or_check_if_view( data, tqa::make_shape(r0.shape()[0],r0.shape()[1],mymesh.size()) );
    (*this)[0] = r0;
    for (size_t u=1; u<mesh().size(); ++u)  (*this)[u] = rhs(mesh()[u]); 
    fill_tail1(mpl::bool_<has_tail>(), rhs); //tail = rhs (domains::infty());  
   }

   TRIQS_NVL_HAS_AUTO_ASSIGN(); 
   template<typename F> friend void triqs_nvl_auto_assign (gf_impl & x, F f) { // mesh is invariant in this case...
    for (size_t u=0; u<x.mesh().size(); ++u)  { x[u] = f(x.mesh()[u]); }
    fill_tail(mpl::bool_<gf_impl::has_tail>(), x.tail, f); 
    // if f is an expression, replace the placeholder with a simple tail. If f is a function callable on infty, 
    // it uses the fact that tail_non_view_type can be caster into domains::infty (See below).
   }

   private : // impl detail. We do not want that f or rhs is called when there is no tail...
   // using this overload, we guarantee that the call of f/rhs will not be *compiled* when not needed... 
   template<typename RHS> static void fill_tail ( mpl::false_, tail_type & t, RHS const & rhs) {} 
   template<typename RHS> static void fill_tail ( mpl::true_,  tail_type & t, RHS const & rhs) {t = rhs( tail::omega(t.shape(),t.size()));} 
   template<typename RHS> void fill_tail1 ( mpl::false_, RHS const & rhs) {} 
   template<typename RHS> void fill_tail1 ( mpl::true_,  RHS const & rhs) { this->tail = rhs (domains::infty()); }

   public : // call operators
   typedef arrays::matrix_view<value_type,       storage_order>  mv_type;
   typedef arrays::matrix_view<const value_type, storage_order>  const_mv_type;

   mv_type       operator() (meshes::mesh_pt<mesh_type> const & x)       { return data(range(),range(),x.i);}
   const_mv_type operator() (meshes::mesh_pt<mesh_type> const & x) const { return data(range(),range(),x.i);}

   tail_view_type       operator() ( domains::infty const & x)       {return tail;}
   const tail_view_type operator() ( domains::infty const & x) const {return tail;}

   typedef typename domain_type::point_type arg0_type;
   mv_type       operator() (arg0_type const & i)       { return this->mesh().interpolate(*this,i);}
   const_mv_type operator() (arg0_type const & i) const { return this->mesh().interpolate(*this,i);}

   /// Save the Green function in i omega_n (as 2 columns).
   void save(std::string file,  bool accumulate=false) const {}

   /// Load the GF
   void load(std::string file){}

   /// 
   friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, gf_impl const & g) {
    BOOST_AUTO( gr , fg.create_group(subgroup_name) );
    h5_write(gr,"data",g.data);
    h5_write(gr,"tail",g.tail);
    //h5_write(gr,"mesh",g.mymesh);
   }

   friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, gf_impl & g){
    BOOST_AUTO( gr,  fg.open_group(subgroup_name) );
    h5_read(gr,"data",g.data);
    h5_read(gr,"tail",g.tail);
    //h5_read(gr,"mesh",g.mymesh);
   }

   //  BOOST Serialization
   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
     ar & boost::serialization::make_nvp("data",data);
     ar & boost::serialization::make_nvp("tail",tail);
     ar & boost::serialization::make_nvp("mesh",mymesh);
    }

   /// print
   friend std::ostream & operator << (std::ostream & out, gf_impl const & x) { return out<<"gf_view";}
 };

 // ---------------------------------------------------------------------------------
 ///The regular class of GF
 template<typename MeshType> class gf :  public gf_impl<MeshType,false> {
  typedef gf_impl<MeshType,false> B;
  public : 
  gf():B() {} 
  gf(size_t N1, size_t N2, MeshType const & m, int tail_order_max = 3):B(N1,N2,m,tail_order_max) {}

  gf(gf const & g): B(g){}
  gf(gf_view<MeshType> const & g): B(g){} 
  template<typename GfType> gf(GfType const & x): B() { *this = x;} 

  using B::operator=;// or the default is = synthetized...

  typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
  typedef gf_view<MeshType> view_type;
  typedef gf<MeshType>      non_view_type;

  using B::operator();
  TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

  /// Slice in orbital space
  const view_type slice(range R1, range R2) const { return view_type(*this,R1,R2); } 

  /// Slice of the mesh
  view_type       slice_mesh(typename MeshType::slice_arg_type arg)       { return view_type (*this, arg);}
  view_type const slice_mesh(typename MeshType::slice_arg_type arg) const { return view_type (*this, arg);}
 };

 // ---------------------------------------------------------------------------------
 ///The View class of GF
 template<typename MeshType> class gf_view :  public gf_impl<MeshType,true> {
  typedef gf_impl<MeshType,true> B;
  public :

  gf_view(gf_view const & g): B(g){}
  gf_view(gf<MeshType> const & g): B(g){}

  gf_view (gf_view const & m,      range R1, range R2) : B(m,R1,R2){} // slice constructor 
  gf_view (gf<MeshType> const & m, range R1, range R2) : B(m,R1,R2){} // slice constructor 

  using B::operator=;// or the default is = synthetized...

  friend std::ostream & triqs_nvl_formal_print(std::ostream & out, gf_view const & x) { return out<<"gf_view";}

  typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 
  typedef gf_view<MeshType> view_type;
  typedef gf<MeshType>      non_view_type;

  using B::operator();
  TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

  /// Slice in orbital space
  const view_type slice(range R1, range R2) const { return view_type(*this,R1,R2); } 

  /// Slice of the mesh : to be implemented
  view_type       slice_mesh(typename MeshType::slice_arg_type arg)       { return view_type (*this, arg);}
  view_type const slice_mesh(typename MeshType::slice_arg_type arg) const { return view_type (*this, arg);}

 };

 // -------------------------------   Expression template for gf  --------------------------------------------------

 // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
 //template <typename T> struct is_scalar_or_element : mpl::or_< tqa::ImmutableMatrix<T>, tup::is_in_ZRC<T> > {};
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

 // grammar and domain and proto expression
 template <typename Expr> struct gf_expr;
 typedef tup::domain_with_copy<gf_grammar,gf_expr> gf_expr_domain;

 template<typename Expr> struct gf_expr : proto::extends<Expr, gf_expr<Expr>, gf_expr_domain> {
  gf_expr( Expr const & expr = Expr() ) : proto::extends<Expr, gf_expr<Expr>, gf_expr_domain> (expr) {}
  typedef typename boost::result_of<gf_mesh_t(Expr,gf_no_mesh) >::type mesh_type;
  mesh_type mesh() const {return gf_mesh_t()(*this,gf_no_mesh());} 
  template<typename T> typename boost::result_of<gf_eval_t(Expr,bf::vector<T>)>::type 
   operator()(T const & x) const {return gf_eval_t()(*this,bf::make_vector(x));}
  friend std::ostream &operator <<(std::ostream &sout, gf_expr<Expr> const &expr) { return tup::print_algebra(sout,expr);} 
 };

 BOOST_PROTO_DEFINE_OPERATORS(LocalGf, gf_expr_domain);
}}}
#endif
