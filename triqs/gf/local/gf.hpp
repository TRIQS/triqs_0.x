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

// PUT AT THE TOP !!!!!
#define BOOST_RESULT_OF_USE_DECLTYPE
#include <boost/utility/result_of.hpp>

#include <boost/type_traits/is_complex.hpp> 
#include <triqs/lazy/core.hpp>
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/utility/proto/tools.hpp>
#include <triqs/arrays/h5/simple_read_write.hpp>
//#include <triqs/utility/proto/algebra.hpp>
#include <triqs/arrays/expressions/matrix_algebra.hpp>
#include <triqs/arrays/expressions/array_algebra.hpp>
#include "./domains.hpp"
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/cat.hpp>

// move it up in the array lib and document it  
namespace triqs { namespace arrays { 

 template< typename A> 
  typename boost::enable_if<is_view_class<A> >::type 
  resize_or_check_if_view ( A & a, typename A::shape_type const & sha) { 
   if (a.shape()!=sha) TRIQS_RUNTIME_ERROR<< "Size mismatch : view class shape = "<<a.shape() << " expected "<<sha;
  }

 template< typename A> 
  typename boost::enable_if<is_value_class<A> >::type 
  resize_or_check_if_view ( A & a, typename A::shape_type const & sha) { if (a.shape()!=sha) a.resize(sha); }

}}

namespace triqs { 

 template<typename MeshType> struct mesh_pt {
  MeshType const & m;
  typename MeshType::index_type i;
  mesh_pt( MeshType const & mesh, typename MeshType::index_type const & index): m(mesh), i(index) {}
  typedef typename MeshType::domain_type::element_type cast_type;
  operator cast_type() const; // cast into the element type of the domain (e.g. real time, real frequency).
 };

 template<typename MeshType> 
  mesh_pt<MeshType> make_mesh_pt(MeshType const & m, typename MeshType::index_type const & i){ return mesh_pt<MeshType>(m,i);}
}

// --------------------------------------------------------------

namespace triqs { namespace gf { namespace local {

 namespace tqa= triqs::arrays; namespace tql= triqs::lazy; namespace mpl= boost::mpl; using tqa::range;
 namespace bf = boost::fusion; namespace tup = triqs::utility::proto; 

 //namespace tag { template<typename MeshType> struct is_gf{};}
 //template <typename G> struct is_gf : boost::is_base_of< tag::is_gf<MeshType>, G> {}; // identification trait

 template<typename MeshType> class gf;
 template<typename MeshType> class gf_view;
 typedef gf     <meshes::tail>      high_frequency_expansion;
 typedef gf_view<meshes::tail>      high_frequency_expansion_view;
 template <typename G> struct is_gf:mpl::false_{}; // trait to identify local gf object. Default is false. Overriden later...

 namespace impl {

  //template<typename LHS, typename RHS, typename Enable=void > struct assignment;

  /*------------------------------------------------------------------------------
   * The implementation class for both the view and the regular gf class
   *-----------------------------------------------------------------------------*/
  template<typename MeshType, bool IsView> class gf_impl { //: tag::is_gf<MeshType> {
   friend class gf_impl<MeshType,!IsView>;
   struct no_tail;
   public:
   typedef MeshType mesh_type;
   typedef typename mesh_type::domain_type domain_type;
   typedef gf_view<mesh_type> view_type;
   typedef gf<mesh_type>      non_view_type;
   typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 

   typedef typename domain_type::gf_result_type value_element_type;

   typedef tqa::array <value_element_type,3,arrays::Option::Fortran>             data_non_view_type;
   typedef typename data_non_view_type::view_type                                data_view_type; 
   typedef typename mpl::if_c<IsView, data_view_type, data_non_view_type>::type  data_type;

   typedef typename mpl::if_c<domain_type::has_tail, high_frequency_expansion, no_tail>::type  tail_non_view_type; 
   typedef typename tail_non_view_type::view_type                                              tail_view_type; 
   typedef typename mpl::if_c<IsView, tail_view_type, tail_non_view_type>::type                tail_type;

   typedef tqa::array<std::string,2> indices_type;

   protected:
   mesh_type mymesh;
   data_type data;
   tail_type tail;
   indices_type ind;

   public:
   domain_type const & domain() const { return mymesh.domain();}
   mesh_type const & mesh() const { return mymesh;}

   data_view_type data_view()             { return data;} // redondant with (range()) ....
   const data_view_type data_view() const { return data;}

   tail_view_type tail_view()             { return tail;} // redondant with ()(domain::infty);
   const tail_view_type tail_view() const { return tail;}

   indices_type const & indices() const   { return ind;}

   protected:
   gf_impl() {} // all arrays of zero size (empty)

   gf_impl (size_t N1, size_t N2, mesh_type const & mesh_, indices_type const & indices_) : 
    mymesh(mesh_), data(N1,N2,mesh().size()), tail(N1,N2,mesh().mesh_tail,indices_), ind(indices_){}

   gf_impl (mesh_type const & d, data_view_type const & dat, tail_view_type const & t, indices_type const & ind_) : 
    data(dat), mymesh(d), tail(t), ind(ind_){}

   gf_impl(gf_impl const & x): mymesh(x.mymesh), data(x.data), tail(x.tail), ind(x.indices()){}
   gf_impl(gf_impl<MeshType,!IsView> const & x): mymesh(x.mymesh), data(x.data), tail(x.tail), ind(x.indices()){} // crucial 

   //does not make sense in the view case... The int is there to prevent any use as copy construction...
   template<typename GfType> gf_impl(GfType const & x, mpl::true_ ) : // with tail
    mymesh(x.mesh()), data(x(arg_range_type())), tail(x(domains::infty())), ind(x.indices()){}

   template<typename GfType> gf_impl(GfType const & x, mpl::false_ ) : // without tail
    mymesh(x.mesh()), data(x(arg_range_type())), ind(x.indices()){}

   template<typename RHS> void operator = (RHS const & rhs) { // first the general version  
    static_assert(is_gf<RHS>::value, "The object does not have the local gf concept");
    mymesh = rhs.mesh();
    BOOST_AUTO( sha , rhs(0).shape());
    BOOST_AUTO( data_sh , tqa::make_shape(sha[0],sha[1],mymesh.size()));
    tqa::resize_or_check_if_view( data, data_sh );
    const size_t Nmax = this->data.shape()[2]; for (size_t u=0; u<Nmax; ++u) data(range(),range(),u) = rhs(u);
    //data = rhs(arg_range_type()) ;// buggy because array algebra is not matrix algebra 
    fill_tail(tail,rhs,domains::infty());
   }

   // quicker version for the class and the view 
   void operator = (gf_impl const & rhs) { mymesh = rhs.mesh(); data = rhs.data; tail = rhs.tail;}
   void operator = (gf_impl<MeshType,!IsView> const & rhs) { mymesh = rhs.mesh(); data = rhs.data; tail = rhs.tail;}

   private: // fill_tail does the calculation only when it makes sense .... 
   template<typename Arg, typename RHS> void fill_tail (high_frequency_expansion &t, RHS const & rhs, Arg const & args ) {t = rhs (args);}
   template<typename Arg, typename RHS> void fill_tail (high_frequency_expansion_view &t, RHS const & rhs, Arg const & args ) {t = rhs (args);}
   template<typename Arg, typename RHS> void fill_tail (no_tail &t, RHS const & rhs, Arg const & args ) {}

   public:
   template<typename F> void set_from_function(F f) { // mesh is invariant in this case... 
    const size_t Nmax = this->data.shape()[2]; for (size_t u=0; u<Nmax; ++u) data(range(),range(),u) = f(1.0*u);
    fill_tail(tail, f ,tail_non_view_type(data.shape()[0], data.shape()[1],mesh().mesh_tail ,indices()));
   }

   // useful ? In the concept ???
   //triqs::arrays::mini_vector<size_t,2> result_shape() const { return data(range(),range(),0).shape();}

   typedef typename mesh_type::index_type arg0_type;
   typedef mesh_pt<mesh_type> mesh_pt_type;
   typedef typename mesh_type::range_type arg_range_type;

   typedef tqa::matrix_view<value_element_type, arrays::Option::Fortran>       mv_type;
   typedef tqa::matrix_view<const value_element_type, arrays::Option::Fortran> const_mv_type;

   template<typename X> struct result;  // implementing ResultOf concept
   template<typename THIS> struct result<THIS(arg0_type)>      { typedef mv_type type; };
   template<typename THIS> struct result<THIS(arg_range_type)> { typedef data_view_type type; };
   template<typename THIS> struct result<THIS(domains::infty)> { typedef tail_view_type type; };

   mv_type       operator() ( mesh_pt_type const & x)       { return data(range(),range(),x.i);}
   const_mv_type operator() ( mesh_pt_type const & x) const { return data(range(),range(),x.i);}
   
   mv_type       operator() ( arg0_type const & i)       { return data(range(),range(),i);}
   const_mv_type operator() ( arg0_type const & i) const { return data(range(),range(),i);}

   data_view_type       operator() ( arg_range_type const & R)       { return data(range(),range(),R);}
   const data_view_type operator() ( arg_range_type const & R) const { return data(range(),range(),R);}

   tail_view_type       operator() ( domains::infty const & x)       { return tail;}
   const tail_view_type operator() ( domains::infty const & x) const { return tail;}

   TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

   template<typename A, typename B> view_type slice(A a, B b) {
    range ra = range(a), rb = range(b); 
    return view_type( gf_impl<mesh_type,true>(mymesh, data (ra,rb, range()), tail.slice(a,b), indices()(ra,rb) )); 
   }
   template<typename A, typename B> const view_type slice(A a, B b) const { 
    range ra = range(a), rb = range(b); 
    return view_type( gf_impl<mesh_type,true>(mymesh, data (ra,rb, range()), tail.slice(a,b), indices()(ra,rb) ));      
   }

   /// Save the Green function in i omega_n (as 2 columns).
   void save(std::string file,  bool accumulate=false) const {}

   /// Load the GF
   void load(std::string file){}

   /// HDF5 saving .... // shall we use boost::file_system for paths ???
   // NOT OK on gcc since CommonFG is abstract... Grrr....
   //friend void h5_write( H5::CommonFG file, std::string path, gf_impl const & g) {} 
   //friend void h5_read( H5::CommonFG file, std::string path, gf_impl & g) {} // exceptions will propagate from the array read/write... 

   private : 
   struct no_tail { // used in place of tail when there is none... an object that does .. no_tail 
    typedef no_tail view_type;
    typedef no_tail non_view_type;
    template<typename A, typename B> view_type slice(A a, B b) const { return no_tail();}
    template<typename T1, typename T2, typename T3, typename T4> no_tail(T1,T2,T3,T4){}
    template<typename T1> no_tail(T1){}
    no_tail(...){}
    template<typename T> void operator=(T) {}
    template<typename T> void set_from_function(T) {}
    void save(std::string file,  bool accumulate=false) const {}
    void load(std::string file){}
    //friend void h5_write( H5::CommonFG file, std::string path, no_tail const & g) {} 
    //friend void h5_read( H5::CommonFG file, std::string path, no_tail & g) {} // exceptions will propagate from the array read/write... 
   };

  };

  // derive from this to add a cast operator to infty only for tail... 
  template<typename MeshType> struct cast_to_infty  {}; 
  template<> struct cast_to_infty<meshes::tail>  { operator domains::infty() const { return domains::infty();} };

  }// namespace impl 

  // --------------------------- THE USER CLASSES ------------------------------------------------------
  ///The regular class of GF
  template<typename MeshType> class gf : public impl::gf_impl<MeshType,false>, public impl::cast_to_infty<MeshType> {
   typedef impl::gf_impl<MeshType,false> B;
   public : 
   gf():B() {} 
   gf (size_t N1, size_t N2, typename B::mesh_type const & m, typename B::indices_type const & i):B(N1,N2,m,i) {}
   template<typename GfType> gf(GfType const & x): B(x,mpl::bool_<B::domain_type::has_tail>()) {}
   template<typename RHS> gf & operator = (RHS const & rhs) { B::operator = (rhs); return *this; } 
  };

  // ---------------------------------------------------------------------------------
  ///The View class of GF
  template<typename MeshType> class gf_view : public impl::gf_impl<MeshType,true>, public impl::cast_to_infty<MeshType> {
   typedef impl::gf_impl<MeshType,true> B;
   public :
   gf_view(gf_view const & g): B(g){}
   gf_view(gf<MeshType> const & g): B(g){}
   // to be removed after changing to auto_assign
   template<typename F> void set_from_function(F f) { B::set_from_function(f);} // bug of autodetection in triqs::lazy on gcc   
   template<typename RHS> gf_view & operator = (RHS const & rhs) { B::operator = (rhs); return *this; } 
   std::ostream & print_for_lazy(std::ostream & out) const { return out<<"gf_view";}
  };

  template <typename M> struct is_gf< gf<M> >:mpl::true_{};
  template <typename M> struct is_gf< gf_view<M> >:mpl::true_{};

  // -------------------------------   Expression template for gf  --------------------------------------------------

  namespace proto=boost::proto;

  // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
  template <typename T> struct is_scalar_or_element : mpl::or_< tqa::expressions::matrix_algebra::IsMatrix<T>, triqs::utility::proto::is_in_ZRC<T> > {};

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
   struct get_mesh {
    BOOST_PROTO_CALLABLE();
    template<typename Sig> struct result;
    template<typename This, typename X> struct result<This(X)> { typedef typename boost::remove_reference<X>::type::mesh_type type;};
    template<typename X> typename X::mesh_type  operator ()(X const & x) const { return x.mesh();}
   };

   struct combine_mesh {
    BOOST_PROTO_CALLABLE();
    template<typename Sig> struct result;
    template<typename This, typename TAG, typename M1, typename M2> struct result<This(TAG,M1,M2)> {typedef M2 type;};
    template<typename This, typename TAG, typename M1> struct result<This(TAG,M1,no_mesh)> {typedef M1 type;};
    template<typename M, typename TAG> M operator ()(TAG, no_mesh const & m1, M const & m2) const { return m2;}
    template<typename M, typename TAG> M operator ()(TAG, M const & m1, no_mesh const & m2) const { return m1;}
    template<typename M1, typename M2, typename TAG> M1 operator ()(TAG, M1 const & m1, M2 const & m2) const { 
     static_assert(boost::is_same<M1,M2>::value, "FATAL : two meshes of different type mixed in an expression");
     return m1;
    }
   };

   struct dom_t : 
    proto::or_<
    proto::when< ScalarGrammar, no_mesh() >
    ,proto::when< ElementGrammar, get_mesh(proto::_value) >
    ,proto::when< proto::plus <dom_t,dom_t>,       combine_mesh (proto::tag::plus(), dom_t(proto::_left), dom_t( proto::_right)) >
    ,proto::when< proto::minus <dom_t,dom_t>,      combine_mesh (proto::tag::plus(), dom_t(proto::_left), dom_t( proto::_right)) >
    ,proto::when< proto::multiplies <dom_t,dom_t>, combine_mesh (proto::tag::plus(), dom_t(proto::_left), dom_t( proto::_right)) >
    ,proto::when< proto::divides <dom_t,dom_t>,    combine_mesh (proto::tag::plus(), dom_t(proto::_left), dom_t( proto::_right)) >
    ,proto::when< proto::unary_expr<proto::_,dom_t >,  dom_t(proto::_left) >
    > {};

   // grammar and domain and proto expression
   template <typename Expr> struct gf_expr;
   typedef tup::domain<gf_grammar,gf_expr,true> gf_expr_domain;

   typedef eval_t eval1;
   template<typename Expr> struct gf_expr : boost::proto::extends<Expr, gf_expr<Expr>, gf_expr_domain>{
    gf_expr( Expr const & expr = Expr() ) : boost::proto::extends<Expr, gf_expr<Expr>, gf_expr_domain>  ( expr ) {}
    typedef typename boost::result_of<dom_t(Expr) >::type mesh_type;
    typedef typename mesh_type::domain_type domain_type;
    mesh_type mesh() const { return dom_t()(*this); } 
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

  // -------------------------------   Expression template for high_frequency_expansion and view -----------------------
#ifdef TAIL_EXPR
 // not written
 //
  namespace tail_expr_templ { 

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

   typedef tqa::matrix<std::complex<double> > matrix_view_type;

   struct tail_mult { 
    BOOST_PROTO_CALLABLE();
    template<typename Sig> struct result;
    template<typename This, typename L, typename R, typename AL> struct result<This(L,R,AL)> { typedef matrix_view_type type;}
    template<typename R, typename L, typename AL>
     matrix_view_type operator ()(L const &l, R const & r, AL const & al) const {
      int n = bf::at<mpl::int_<0> >(al);//get n
      //if (n1!=t.n1 || n2!=t.n2) TRIQS_RUNTIME_ERROR<<"multiplication is valid only for similar tail shapes !";
      //if (new_ordermin < OrderMinMIN) TRIQS_RUNTIME_ERROR<<"The multiplication makes the new tail have a too small OrderMin";
      matrix_view_type::non_view_type res(l(l.mesh().order_min()).shape()[0], l(l.mesh().order_min()).shape()[1]); res()=0;
      const int kmin = std::max(0, n - r.mesh().order_max() - l.mesh().order_min() );
      const int kmax = std::min(l.mesh().order_max() - l.mesh().order_min(), n - r.mesh().order_min() - l.mesh().order_
	min() );
      std::cout<< kmin << "  "<<kmax<<std::endl;
      for (int k = kmin; k <= kmax; ++k)  { std::cout <<" k = "<< res<<std::endl; 
       res += l(l.mesh().order_min() +k) * r( n- l.mesh().order_min() -k);}
       return res;
     }
   };

   struct eval_t : 
    proto::or_<
    proto::when<ScalarGrammar, proto::_value>
    ,proto::when<ElementGrammar, tup::eval_fnt<1>(proto::_value, proto::_state) >
    ,proto::when< proto::plus <eval_t,eval_t>,       proto::_default<eval_t > >
    ,proto::when< proto::minus <eval_t,eval_t>,      proto::_default<eval_t > >
    ,proto::when< proto::multiplies <eval_t,eval_t>, tail_mult(proto::_left, proto::_right, proto::state) >
    ,proto::when< proto::divides <eval_t,eval_t>,    proto::_default<eval_t > >
    ,proto::when< proto::unary_expr<proto::_,eval_t >,  eval_t(proto::_left) >
    ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::_default<eval_t > >
    > {};

   // computation of domain ....
   template<typename L, typename R> struct multiplies { 
    L const & l; R const & r; multiplies (L const & l_, R const & r_):l(l_),r(r_) {}
    typedef meshes::tail  mesh_type;
    meshes::tail domain () const { 
     int omin = l.order_min()+r.order_min(); 
     return meshes::tail(omin,omin + std::min(l.domain().len(),r.mesh().len()));
    }

    struct no_mesh{ typedef void domain_type;}; 
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
      static_assert(boost::is_same<M1,M2>::value, "FATAL : two meshes of different type mixed in an expression");
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

    // grammar and domain and proto expression
    template <typename Expr> struct tail_expr;
    typedef tup::domain<gf_grammar,tail_expr,true> tail_expr_domain;

    typedef eval_t eval1;
    template<typename Expr> struct tail_expr : boost::proto::extends<Expr, tail_expr<Expr>, tail_expr_domain>{
     tail_expr( Expr const & expr = Expr() ) : boost::proto::extends<Expr, tail_expr<Expr>, tail_expr_domain>  ( expr ) {}
     typedef typename boost::result_of<dom_t(Expr) >::type mesh_type;
     typedef typename mesh_type::domain_type domain_type;
     mesh_type mesh() const { return dom_t()(*this); } 
     template<typename T> typename boost::result_of<eval1(Expr,bf::vector<T>) >::type 
      operator() (T const & x) const { return eval1()(*this, bf::make_vector(x)); }
     typedef tqa::array_view<std::string,2> indices_type;
     indices_type indices() const { return indices_type::non_view_type();}
     friend std::ostream &operator <<(std::ostream &sout, tail_expr<Expr> const &expr) { return boost::proto::eval(expr, tup::algebra_print_ctx (sout)); }
    };
    template <typename Expr> struct is_gf<tail_expr< Expr > >:mpl::true_{};

    BOOST_PROTO_DEFINE_OPERATORS(is_tail, tail_expr_domain);
   }
#endif

  }}}

#endif

