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

namespace triqs { namespace gf { namespace local {

 namespace tqa= triqs::arrays; namespace tql= triqs::lazy; namespace mpl= boost::mpl; using tqa::range;
 namespace bf = boost::fusion; namespace tup = triqs::utility::proto; 

 //namespace tag { template<typename DomainType> struct is_gf{};}
 //template <typename G> struct is_gf : boost::is_base_of< tag::is_gf<DomainType>, G> {}; // identification trait

 template<typename DomainType> class gf;
 template<typename DomainType> class gf_view;
 typedef gf<meshes::tail>      high_frequency_expansion;
 typedef gf_view<meshes::tail> high_frequency_expansion_view;
 template <typename G> struct is_gf:mpl::false_{}; // trait to identify local gf object. Default is false. Overriden later...

 namespace impl {

  //template<typename LHS, typename RHS, typename Enable=void > struct assignment;

  /*------------------------------------------------------------------------------
   * The implementation class for both the view and the regular gf class
   *-----------------------------------------------------------------------------*/
  template<typename DomainType, bool IsView> class gf_impl { //: tag::is_gf<DomainType> {
   friend class gf_impl<DomainType,!IsView>;
   struct no_tail;
   public:
   typedef DomainType domain_type;
   typedef gf_view<domain_type> view_type;
   typedef gf<domain_type>      non_view_type;
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
   domain_type dom;
   data_type data;
   tail_type tail;
   indices_type ind;

   gf_impl() {} // all arrays of zero size (empty)

   gf_impl (size_t N1, size_t N2, domain_type const & domain_, indices_type const & indices_) : 
    dom(domain_), data(N1,N2,domain().size()), tail(N1,N2,domain().mesh_tail,indices_), ind(indices_){}

   gf_impl (domain_type const & d, data_view_type const & dat, tail_view_type const & t, indices_type const & ind_) : 
    data(dat), dom(d), tail(t), ind(ind_){}

   //gf_impl(gf_impl const & x): dom(x.dom), data(x.data), tail(x.tail), ind(x.indices()){}

   //template<typename GfType> gf_impl(GfType const & x): dom(x.domain()), data(x.data_view()), tail(x.tail_view()), ind(x.indices()){}
   //template<typename GfType> gf_impl(GfType const & x): dom(x.domain()), data(x(arg_range_type())), tail(x(domains::infty())), ind(x.indices()){}
   template<typename GfType> gf_impl(GfType const & x): dom(x.domain()), data(x(arg_range_type())), tail(call_at_infty<tail_view_type>(x)), ind(x.indices()){}

   template<typename RHS> void operator = (RHS const & rhs) { 
    static_assert(is_gf<RHS>::value, "The object does not have the local gf concept");
    dom = rhs.domain(); data = rhs(arg_range_type()) ; tail= call_at_infty<tail_view_type>(rhs); //(domains::infty()); 
   } 

   private : 
   template<typename TailViewType, typename F> TailViewType call_at_infty(F const & f) { return typename TailViewType::non_view_type( f(domains::infty()));}
   template<typename TailViewType, typename F> TailViewType call_at_infty(F const & f, TailViewType const & arg) { return typename TailViewType::non_view_type( f(arg));}
   template<typename F> no_tail call_at_infty(F const & f) { return no_tail();}
   template<typename F> no_tail call_at_infty(F const & f,no_tail const &) { return no_tail();}

   public:

   domain_type const & domain() const { return dom;}

   data_view_type data_view()             { return data;} // redondant with (range()) ....
   const data_view_type data_view() const { return data;}

   tail_view_type tail_view()             { return tail;} // redondant with ()(domain::infty);
   const tail_view_type tail_view() const { return tail;}

   indices_type const & indices() const   { return ind;}

   // useful ? In the concept ???
   //triqs::arrays::mini_vector<size_t,2> result_shape() const { return data.shape();}

   typedef typename domain_type::index_type arg0_type;
   typedef typename domain_type::range_type arg_range_type;

   typedef tqa::matrix_view<value_element_type, arrays::Option::Fortran>       mv_type;
   typedef tqa::matrix_view<const value_element_type, arrays::Option::Fortran> const_mv_type;

   template<typename X> struct result;                      // implementing ResultOf concept
   template<typename THIS> struct result<THIS(arg0_type)>      { typedef mv_type type; };
   template<typename THIS> struct result<THIS(arg_range_type)> { typedef typedef data_view_type type; };
   template<typename THIS> struct result<THIS(domains::infty)> { typedef typedef tail_view_type type; };

   mv_type       operator() ( arg0_type const & x)       { return data(range(),range(),x);}
   const_mv_type operator() ( arg0_type const & x) const { return data(range(),range(),x);}

   data_view_type       operator() ( arg_range_type const & R)       { return data(range(),range(),R);}
   const data_view_type operator() ( arg_range_type const & R) const { return data(range(),range(),R);}

   tail_view_type       operator() ( domains::infty const & x)       { return tail;}
   const tail_view_type operator() ( domains::infty const & x) const { return tail;}

   TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

   template<typename A, typename B> view_type slice(A a, B b) {
    range ra = range(a), rb = range(b); 
    return view_type( gf_impl<domain_type,true>(dom, data (ra,rb, range()), tail.slice(a,b), indices()(ra,rb) )); 
   }
   template<typename A, typename B> const view_type slice(A a, B b) const { 
    range ra = range(a), rb = range(b); 
    return view_type( gf_impl<domain_type,true>(dom, data (ra,rb, range()), tail.slice(a,b), indices()(ra,rb) ));      
   }

   template<typename F> void set_from_function(F f) { // domain is invariant in this case... 
    const size_t Nmax = this->data.shape()[2]; for (size_t u=0; u<Nmax; ++u) this->data(range(),range(),u) = f(1.0*u);
    //tail = call_at_infty(f ,tail_non_view_type(data.shape()[0], data.shape()[1],domain().mesh_tail ,indices()));
    //compute_tail(tail,f); // the tail is only computed if it exists, or f maybe called for infty in undefined cases.
   }
   private:
   template<typename F> void compute_tail(gf_impl<meshes::tail, IsView> & t, F f) {
    // if f is only defined for domain::infty, the cast operator added in tail will do the job... 
    t = call_at_infty(f ,tail_non_view_type(data.shape()[0], data.shape()[1],domain().mesh_tail ,indices()));
    //t = f( tail_non_view_type(data.shape()[0], data.shape()[1],domain().mesh_tail ,indices()));
    //t = f( wrap_infty<tail_view_type,domain_type>( tail_non_view_type(data.shape()[0], data.shape()[1],domain().mesh_tail ,indices())));
   }
   template<typename F> void compute_tail(no_tail & t, F f) {}

   public:
   /// Save the Green function in i omega_n (as 2 columns).
   void save(std::string file,  bool accumulate=false) const {}

   /// Load the GF
   void load(std::string file){}

   /// HDF5 saving .... // shall we use boost::file_system for paths ???
   friend void h5_write( H5::CommonFG file, std::string path, gf_impl const & g) {} 
   friend void h5_read( H5::CommonFG file, std::string path, gf_impl & g) {} // exceptions will propagate from the array read/write... 

   private : 
   struct no_tail { // used in place of tail when there is none... an object that does .. no_tail 
    typedef no_tail view_type;
    typedef no_tail non_view_type;
    template<typename A, typename B> view_type slice(A a, B b) const { return no_tail();}
    no_tail(...){}
    template<typename T1, typename T2, typename T3, typename T4> no_tail(T1,T2,T3,T4){}
    template<typename T> void operator=(T) {}
    template<typename T> void set_from_function(T) {}
   };

  };

  // derive from this to add a cast operator to infty only for tail... 
  template<typename DomainType> struct cast_to_infty  {}; 
  template<> struct cast_to_infty<meshes::tail>  { operator domains::infty() const { return domains::infty();} };

  }// namespace impl 

  // --------------------------- THE USER CLASSES ------------------------------------------------------
  ///The regular class of GF
  template<typename DomainType> class gf : public impl::gf_impl<DomainType,false>, public impl::cast_to_infty<DomainType> {
   typedef impl::gf_impl<DomainType,false> base_type;
   public : 
   gf():base_type() {} 
   gf (size_t N1, size_t N2, typename base_type::domain_type const & domain_, typename base_type::indices_type const & indices_) : 
    base_type(N1,N2,domain_,indices_) {}
   template<typename GfType> gf(GfType const & x): base_type(x) {}
   template<typename RHS> gf & operator = (RHS const & rhs) { base_type::operator = (rhs); return *this; } 
  };

  // ---------------------------------------------------------------------------------
  ///The View class of GF
  template<typename DomainType> class gf_view : public impl::gf_impl<DomainType,true>, public impl::cast_to_infty<DomainType> {
   typedef impl::gf_impl<DomainType,true> base_type;
   public :
   template<typename GfType> gf_view(GfType const & x): base_type(x) {};
   // to be removed after changing to auto_assign
   template<typename F> void set_from_function(F f) { base_type::set_from_function(f);} // bug of autodetection in triqs::lazy on gcc   
   template<typename RHS> gf_view & operator = (RHS const & rhs) { base_type::operator = (rhs); return *this; } 
   std::ostream & print_for_lazy(std::ostream & out) const { return out<<"gf_view";}
  };

  template <typename M> struct is_gf< gf<M> >:mpl::true_{};
  template <typename M> struct is_gf< gf_view<M> >:mpl::true_{};

  // -------------------------------   Expression template for each domain --------------------------------------------------

  namespace proto=boost::proto;

  // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
  template <typename T> struct is_scalar_or_element : mpl::or_< tqa::expressions::matrix_algebra::IsMatrix<T>, triqs::utility::proto::is_in_ZRC<T> > {};


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


  struct eval_transform : 
   proto::or_<
   proto::when<ScalarGrammar, proto::_value>
   ,proto::when<ElementGrammar, tup::eval_fnt<1>(proto::_value, proto::_state) >
   ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::_default<eval_transform > >
   > {};

  /*
 struct node_desc_impl_tail : node_desc_impl <meshes::tail> {

  template<typename L, typename R> struct multiplies { 
   L const & l; R const & r; multiplies (L const & l_, R const & r_):l(l_),r(r_) {}
   typedef meshes::tail  domain_type;
   meshes::tail domain () const { 
    int omin = l.order_min()+r.order_min(); 
    return meshes::tail(omin,omin + std::min(l.domain().len(),r.domain().len()));
   }

   tail_result_type operator() (int n) const {
    //if (n1!=t.n1 || n2!=t.n2) triqs_runtime_error<<"multiplication is valid only for similar tail shapes !";
    //if (new_ordermin < OrderMinMIN) TRIQS_RUNTIME_ERROR<<"The multiplication makes the new tail have a too small OrderMin";
    tail_result_type::non_view_type res(l(l.domain().order_min()).shape()[0], l(l.domain().order_min()).shape()[1]); res()=0;
    const int kmin = std::max(0, n - r.domain().order_max() - l.domain().order_min() );
    const int kmax = std::min(l.domain().order_max() - l.domain().order_min(), n - r.domain().order_min() - l.domain().order_
min() );
    std::cout<< kmin << "  "<<kmax<<std::endl;
    for (int k = kmin; k <= kmax; ++k)  { std::cout <<" k = "<< res<<std::endl; res += l(l.domain().order_min() +k) * r( n- l
.domain().order_min() -k);}
    return res;
   }
  };
 };
*/

  struct no_mesh{}; 
  //template<typename A, typename Enable=void> struct get_mesh_type { typedef no_mesh type;};
  //template<typename A> struct get_mesh_type<A,typename boost::enable_if<is_gf<A> >::type > { typedef typename A::domain_type;};

  /*  struct get_mesh_scalar {
      BOOST_PROTO_CALLABLE();
      template<typename Sig> struct result;
      template<typename This, typename X> struct result<This(X)> { typedef no_mesh type;};
      template<typename X> typename operator ()(X) const { return no_mesh();}
      };
      */
  struct get_mesh {
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename X> struct result<This(X)> { typedef typename boost::remove_reference<X>::type::domain_type type;};
   template<typename X> typename X::domain_type  operator ()(X const & x) const { return x.domain();}
  };

  struct combine_mesh {
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename M1, typename M2> struct result<This(M1,M2)> {typedef M2 type;};
   template<typename This, typename M1> struct result<This(M1,no_mesh)> {typedef M1 type;};
   //template<typename This, typename M1, typename M2> struct result<This(M1,M2)> : mpl::if_<boost::is_same<M1,no_mesh>, M2, M1> {};
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
  template <typename Expr> struct gf_expr;
  typedef tup::domain<gf_grammar,gf_expr,true> gf_expr_domain;

  //typedef eval_transform<is_gf, is_scalar_or_element,eval_fnt<1> > eval1;
  typedef eval_transform eval1;
  template<typename Expr> struct gf_expr : boost::proto::extends<Expr, gf_expr<Expr>, gf_expr_domain>{
   gf_expr( Expr const & expr = Expr() ) : boost::proto::extends<Expr, gf_expr<Expr>, gf_expr_domain>  ( expr ) {}
   typedef typename boost::result_of<dom_t(Expr) >::type domain_type;
   domain_type domain() const { return dom_t()(*this); } //, no_mesh()); }
   template<typename T> 
    typename boost::result_of<eval1(Expr,bf::vector<T>) >::type operator() (T const & x) const {return eval1()(*this, bf::make_vector(x));}
   typedef tqa::array_view<std::string,2> indices_type;
   indices_type indices() const { return indices_type::non_view_type();}
   friend std::ostream &operator <<(std::ostream &sout, gf_expr<Expr> const &expr) { return boost::proto::eval(expr, tup::algebra_print_ctx (sout)); }
 };
 template <typename Expr> struct is_gf<gf_expr< Expr > >:mpl::true_{};

 BOOST_PROTO_DEFINE_OPERATORS(is_gf, gf_expr_domain);

 }}}

#endif

