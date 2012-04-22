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
#include <triqs/utility/proto/algebra.hpp>
#include <triqs/arrays/expressions/matrix_algebra.hpp>
#include "./domains.hpp"
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/cat.hpp>

namespace triqs { namespace gf { namespace local {

 namespace tqa= triqs::arrays; namespace tql= triqs::lazy; namespace mpl= boost::mpl; using tqa::range;
 namespace tup = triqs::utility::proto; 

 namespace tag { template<typename DomainType> struct is_gf{};}

 template<typename DomainType> class gf;
 template<typename DomainType> class gf_view;
 typedef gf<meshes::tail>      high_frequency_expansion;
 typedef gf_view<meshes::tail> high_frequency_expansion_view;

 namespace impl {

  //template<typename LHS, typename RHS, typename Enable=void > struct assignment;

  struct nothing { // used in place of tail when there is none... an object that does .. nothing 
   typedef nothing view_type;
   template<typename A, typename B> view_type slice(A a, B b) const { return nothing();}
   nothing(...){}
   template<typename T1, typename T2, typename T3, typename T4> nothing(T1,T2,T3,T4){}
   template<typename T> void operator=(T) {}
   template<typename T> void set_from_function(T) {}
  };

  template<typename TailViewType, typename DomainType>
   struct wrap_infty  : tag::is_gf<DomainType>  {
    TailViewType tail;
    wrap_infty(TailViewType const & t): tail(t) {}
    operator domains::infty() const { return domains::infty();}
    typename TailViewType::const_result_type operator() ( typename TailViewType::arg0_type const & x) const { return tail(x);}
    typedef meshes::tail domain_type;
   };

  /*------------------------------------------------------------------------------
   * The implementation class for both the view and the regular gf class
   *-----------------------------------------------------------------------------*/
  template<typename DomainType, bool IsView> class gf_impl : tag::is_gf<DomainType> {
    friend class gf_impl<DomainType,!IsView>;
   public:
    typedef DomainType domain_type;
    typedef typename domain_type::gf_result_type value_element_type;
    typedef tqa::array<std::string,2> indices_type;

    typedef gf_view<domain_type> view_type;
    typedef gf<domain_type>      non_view_type;
    typedef void has_view_type_tag; // so that triqs::lazy can put view in the expression trees. 

    typedef tqa::array      <value_element_type,3,arrays::Option::Fortran>        data_non_view_type;
    typedef tqa::array_view <value_element_type,3,arrays::Option::Fortran>        data_view_type;
    typedef typename mpl::if_c<IsView, data_view_type, data_non_view_type>::type  data_type;

    typedef typename mpl::if_c<domain_type::has_tail, high_frequency_expansion, nothing>::type       tail_non_view_type; 
    typedef typename mpl::if_c<domain_type::has_tail, high_frequency_expansion_view, nothing>::type  tail_view_type; 
    typedef typename mpl::if_c<IsView, tail_view_type, tail_non_view_type>::type       tail_type;

   protected:
    domain_type dom;
    data_type data;
    tail_type tail;
    indices_type ind;

   public:

    domain_type const & domain() const { return dom;}

    data_view_type data_view()             { return data;}
    const data_view_type data_view() const { return data;}

    tail_view_type tail_view()             { return tail;}
    const tail_view_type tail_view() const { return tail;}

    indices_type const & indices() const { return ind;}

   protected:
    gf_impl (domain_type const & d, data_view_type const & dat, tail_view_type const & t, indices_type const & ind_) : 
     data(dat), dom(d), tail(t), ind(ind_){}

    gf_impl(gf_impl const & x): dom(x.dom), data(x.data), tail(x.tail), ind(x.indices()){}

    template<typename GfType> gf_impl(GfType const & x): dom(x.domain()), data(x.data_view()), tail(x.tail_view()), ind(x.indices()){}

    gf_impl() {} 

    gf_impl (size_t N1, size_t N2, domain_type const & domain_, indices_type const & indices_) : 
     dom(domain_), data(N1,N2,domain().len()), tail(N1,N2,domain().mesh_tail,indices_), ind(indices_){}

   public:

    triqs::arrays::mini_vector<size_t,2> result_shape() const { return data.shape();}

    typedef typename domain_type::index_type arg0_type;
    typedef tqa::matrix_view<value_element_type, arrays::Option::Fortran>       result_type;
    typedef tqa::matrix_view<const value_element_type, arrays::Option::Fortran> const_result_type;

    result_type       operator() ( arg0_type const & x)       { return data(range(),range(),x);}
    const_result_type operator() ( arg0_type const & x) const { return data(range(),range(),x);}

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

    template<typename RHS> void operator = (RHS const & rhs) { dom = rhs.domain(); data = rhs.data; tail= rhs.data; } 

    template<typename F> void set_from_function(F f) { 
     const size_t Nmax = this->data.shape()[2]; for (size_t u=0; u<Nmax; ++u) this->data(range(),range(),u) = f(1.0*u);
     compute_tail(tail,f); // the tail is only computed if it exists, or f maybe called for infty in undefined cases.
    }
   private:
    template<typename F> void compute_tail(gf_impl<meshes::tail, IsView> & t, F f) { 
     t.set_from_function( f( wrap_infty<tail_view_type,domain_type>( tail_non_view_type(data.shape()[0], data.shape()[1],domain().mesh_tail ,indices()))));
    }
    template<typename F> void compute_tail(nothing & t, F f) {}

   public:
    /// Save the Green function in i omega_n (as 2 columns).
    void save(std::string file,  bool accumulate=false) const {}

    /// Load the GF
    void load(std::string file){}

    /// HDF5 saving ....

  };

 }// namespace impl 

 // ---------------------------------------------------------------------------------
 ///The regular class of GF
 template<typename DomainType> class gf : public impl::gf_impl<DomainType,false> {
  typedef impl::gf_impl<DomainType,false> base_type;
  public : 
  gf():base_type() {} 
  gf (size_t N1, size_t N2, typename base_type::domain_type const & domain_, typename base_type::indices_type const & indices_) : 
   base_type(N1,N2,domain_,indices_) {}
  template<typename GfType> gf(GfType const & x): base_type(x) {}
  template<typename RHS> gf & operator = (RHS const & rhs) { base_type::operator = (rhs); return *this; } 
 };

 // ---------------------------------------------------------------------------------
 /// The View class of GF
 template<typename DomainType> class gf_view : public impl::gf_impl<DomainType,true> {
  typedef impl::gf_impl<DomainType,true> base_type;
  public :
  template<typename GfType> gf_view(GfType const & x): base_type(x) {};
  template<typename F> void set_from_function(F f) { base_type::set_from_function(f);} // bug of autodetection in triqs::lazy on gcc   
  template<typename RHS> gf_view & operator = (RHS const & rhs) { base_type::operator = (rhs); return *this; } 
  std::ostream & print_for_lazy(std::ostream & out) const { return out<<"gf_view";}
 };

 // -------------------------------   Expression template for each domain --------------------------------------------------

 namespace proto=boost::proto;

 // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
 template <typename T> struct is_scalar_or_element : mpl::or_< tqa::expressions::matrix_algebra::IsMatrix<T>, triqs::utility::proto::is_in_ZRC<T> > {};

 template< typename DomainType> struct node_desc_impl {// algebra node description  

  struct no_domain{}; 
  template<typename L, typename R>  static DomainType compute_domain (L const & dl, R const & dr) 
   //{ if (!(dl==dr)) TRIQS_RUNTIME_ERROR <<"domain mismatch"; return dl; }
  {  return dl; }
  template<typename L> static DomainType compute_domain (L const & dl, no_domain) { return dl; }
  template<typename R> static DomainType compute_domain (no_domain, R const & dr) { return dr; }

  template<typename S> struct scalar {
   S s; scalar( S const & x) : s(x) {}
   template<typename T> S operator() (T) const {return s;}
   typedef no_domain domain_type;
   no_domain domain () const {} 
  };
  template<typename L, typename R> struct plus { 
   L const & l; R const & r; plus (L const & l_, R const & r_):l(l_),r(r_) {}
   template<typename T> auto operator() (T const & arg) const -> decltype( l(arg) + r(arg))  {return l(arg)+r(arg);}
   typedef DomainType domain_type;
   DomainType domain () const { return compute_domain(l.domain(),r.domain());}
  }; 
  template<typename L, typename R> struct minus { 
   L const & l; R const & r; minus (L const & l_, R const & r_):l(l_),r(r_) {}
   template<typename T> auto operator() (T const & arg) const -> decltype( l(arg) - r(arg))  {return l(arg)-r(arg);}
   typedef DomainType domain_type;
   DomainType domain () const { return compute_domain(l.domain(),r.domain());}
  }; 
  template<typename L, typename R> struct multiplies { 
   L const & l; R const & r; multiplies (L const & l_, R const & r_):l(l_),r(r_) {}
   template<typename T> auto operator() (T const & arg) const -> decltype( l(arg) * r(arg))  {return l(arg)*r(arg);}
   typedef DomainType domain_type;
   DomainType domain () const { return compute_domain(l.domain(),r.domain());}
  }; 
  template<typename L, typename R> struct divides { 
   L const & l; R const & r; divides (L const & l_, R const & r_):l(l_),r(r_) {}
   template<typename T> auto operator() (T const & arg) const -> decltype( l(arg) / r(arg))  {return l(arg)/r(arg);}
   typedef DomainType domain_type;
   DomainType domain () const { return compute_domain(l.domain(),r.domain());}
  }; 
  template<typename L> struct negate  { 
   L const & l; negate (L const & l_):l(l_) {} 
   template<typename T> auto operator() (T const & arg) const -> decltype( - l(arg) )  {return - l(arg);}
   typedef DomainType domain_type;
   DomainType domain () const { return l.domain(); }
  };
 };

 // a special case for tails
 typedef tqa::matrix_view<std::complex<double> , tqa::Option::Fortran> tail_result_type;

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
    const int kmax = std::min(l.domain().order_max() - l.domain().order_min(), n - r.domain().order_min() - l.domain().order_min() );
    std::cout<< kmin << "  "<<kmax<<std::endl;
    for (int k = kmin; k <= kmax; ++k)  { std::cout <<" k = "<< res<<std::endl; res += l(l.domain().order_min() +k) * r( n- l.domain().order_min() -k);}
    return res;
   }
  };
 };

 // gathering everyone.... 
 template<typename DomainType> struct node_desc : node_desc_impl <DomainType> {};
 template<> struct node_desc<meshes::tail> : node_desc_impl_tail {};

 template< typename DomainType> struct expr_templ {  // now the boost::proto business... 
  template <typename G> struct is_gf : boost::is_base_of< tag::is_gf<DomainType>, G> {}; // identification trait
  // grammar and domain and proto expression
  template <typename Expr> struct gf_expr;

//-------------
//#define OUT_OF_BOX
#ifdef OUT_OF_BOX
  typedef node_desc<DomainType> OpsCompound;
  struct LeafGrammar   : proto::and_< proto::terminal<proto::_>, proto::if_<is_gf<proto::_value>()> > {}; 
    struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<is_scalar_or_element<proto::_value>()> > {}; 

    struct Grammar : 
     proto::or_<
     proto::when< ScalarGrammar,                       typename OpsCompound::template scalar<proto::_value>(proto::_value) >
     ,proto::when< LeafGrammar,                        proto::_value >
     ,proto::when< proto::plus <Grammar,Grammar>,      typename OpsCompound::template plus<proto::_left,proto::_right > (proto::_left,proto::_right) >
     ,proto::when< proto::minus <Grammar,Grammar>,     typename OpsCompound::template minus<proto::_left,proto::_right > (proto::_left,proto::_right)  >
     ,proto::when< proto::multiplies<Grammar,Grammar>, typename OpsCompound::template multiplies<proto::_left,proto::_right > (proto::_left,proto::_right)>
     ,proto::when< proto::divides<Grammar,Grammar>,    typename OpsCompound::template divides<proto::_left,proto::_right > (proto::_left,proto::_right)>
     ,proto::when< proto::negate<Grammar >,            typename OpsCompound::template negate <proto::_left >(proto::_left) >
     > {};


  typedef Grammar grammar;
#else
  typedef typename tup::algebra::grammar_generator< node_desc<DomainType>,is_gf, is_scalar_or_element >::type grammar;
#endif
  typedef tup::domain<grammar,gf_expr,true> gf_domain;

  template<typename Expr> struct gf_expr : boost::proto::extends<Expr, gf_expr<Expr>, gf_domain>{
   gf_expr( Expr const & expr = Expr() ) : boost::proto::extends<Expr, gf_expr<Expr>, gf_domain>  ( expr ) {}
   typedef typename boost::remove_reference<typename boost::result_of<grammar(Expr) >::type>::type _G;
   template<typename T> typename boost::result_of<_G(T)>::type operator() (T const & x) const { return grammar()(*this)(x); }
   typename _G::domain_type domain() const { return grammar()(*this).domain(); }
   friend std::ostream &operator <<(std::ostream &sout, gf_expr<Expr> const &expr) { return boost::proto::eval(expr, tup::algebra::print_ctx (sout)); }
  };
 };

 // BOOST macro for all domains...
#define AUX(r,data,DOM)  BOOST_PROTO_DEFINE_OPERATORS(expr_templ<meshes::DOM>::is_gf, expr_templ<meshes::DOM>::gf_domain);
 BOOST_PP_SEQ_FOR_EACH(AUX, nil , TRIQS_LOCAL_GF_DOMAIN_LIST);
#undef AUX
 BOOST_PROTO_DEFINE_OPERATORS(expr_templ<meshes::tail>::is_gf, expr_templ<meshes::tail>::gf_domain);

}}}

#endif

