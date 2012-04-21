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
 namespace tup = triqs::utility::proto; namespace p_tag= boost::proto::tag;

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
   };

  /*------------------------------------------------------------------------------
   * The implementation class for both the view and the regular gf class
   *-----------------------------------------------------------------------------*/
  template<typename DomainType, bool IsView> class gf_impl : tag::is_gf<DomainType> {

   public:

    typedef DomainType domain_type;
    typedef typename domain_type::gf_result_type value_element_type;
    typedef std::vector<std::vector<std::string> > indices_type;

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
     return view_type( gf_impl<domain_type,true>(dom, data (ra,rb, range()), tail.slice(a,b), slice_indices(indices(),ra,rb) )); 
    }
    template<typename A, typename B> const view_type slice(A a, B b) const { 
     range ra = range(a), rb = range(b); 
     return view_type( gf_impl<domain_type,true>(dom, data (ra,rb, range()), tail.slice(a,b), slice_indices(indices(),ra,rb) )); 
    }

    template<typename RHS> void operator = (RHS const & rhs) { dom = rhs.domain(); data = rhs.data; tail= rhs.data; } 

    // lazy_assignable // TO DO : the computation of the tail is not ready
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

   private : 
    indices_type slice_indices(indices_type const & V, range const & R1, range const & R2){
     indices_type res(2); res[0].reserve(R1.last() - R1.first() + 1); res[1].reserve(R2.last() - R2.first() + 1); 
     for (size_t u =R1.first(); u<=R1.last(); ++u) res[0].push_back(V[0][u]);
     for (size_t u =R2.first(); u<=R2.last(); ++u) res[1].push_back(V[1][u]);
     return res;
    }
    friend class gf_impl<domain_type,!IsView>;
  };

 }// namespace impl 

 // ---------------------------------------------------------------------------------
 /**
  * The View class of GF
  */
 template<typename DomainType> class gf_view : public impl::gf_impl<DomainType,true> {
  typedef impl::gf_impl<DomainType,true> base_type;
  public :

  // gf (domain_type const & domain_, data_view_type const & data_, tail_view_type const & tail_, indices_type const & indices_) : 
  // B(domain_,data_,tail_,indices_) {}

  template<typename GfType> gf_view(GfType const & x): base_type(x) {};

  template<typename F> void set_from_function(F f) { base_type::set_from_function(f);} // bug of autodetection in triqs::lazy on gcc   

  template<typename RHS> gf_view & operator = (RHS const & rhs) { base_type::operator = (rhs); return *this; } 

  std::ostream & print_for_lazy(std::ostream & out) const { return out<<"gf_view";}
 };

 // ---------------------------------------------------------------------------------
 /**
  * The regular class of GF
  */
 template<typename DomainType> class gf : public impl::gf_impl<DomainType,false> {
  typedef impl::gf_impl<DomainType,false> base_type;
  public : 

  gf():base_type() {} 

  gf (size_t N1, size_t N2, typename base_type::domain_type const & domain_, typename base_type::indices_type const & indices_) : 
   base_type(N1,N2,domain_,indices_) {}

  template<typename GfType> gf(GfType const & x): base_type(x) {}

  template<typename RHS> gf & operator = (RHS const & rhs) { base_type::operator = (rhs); return *this; } 

 };

 // -------------------------------   Expression template for each domain --------------------------------------------------

 // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
 template <typename T> struct is_scalar_or_element : mpl::or_< tqa::expressions::matrix_algebra::IsMatrix<T>, triqs::utility::proto::is_in_ZRC<T> > {};
 
#define AUX(r,data,DOM)\
 template <typename G> struct BOOST_PP_CAT(is_gf_,DOM) : boost::is_base_of< tag::is_gf<meshes::DOM>, G> {};\
 TRIQS_PROTO_DEFINE_ALGEBRA_VALUED_FNT_ALG (BOOST_PP_CAT(is_gf_,DOM), is_scalar_or_element);
 
 BOOST_PP_SEQ_FOR_EACH(AUX, nil , TRIQS_LOCAL_GF_DOMAIN_LIST);
#undef AUX

 /***************************************************************************
  *   Computation of domain
  ***************************************************************************/
 template<typename DomainType>
 struct domain_ctx : boost::proto::callable_context< domain_ctx<DomainType> const > {
  typedef DomainType result_type; typedef boost::proto::tag::terminal term_tag;
  template <typename T> typename boost::enable_if <tup::is_in_ZRC<T>, result_type>::type operator ()(term_tag, const T & A) const { return result_type(); }
  template <typename T> typename boost::disable_if<tup::is_in_ZRC<T>, result_type>::type operator ()(term_tag, const T & A) const { return A.domain();}
  template<typename TAG, typename L, typename R> result_type operator ()(TAG, L const &l, R const &r) const { return l.domain(); }
 };

 /***************************************************************************
  *    Implementation of tails 
  ***************************************************************************/
 /*class high_frequency_expansion_view : public impl::gf_impl<meshes::tail,true> {
  typedef impl::gf_impl<meshes::tail,true> base_type;
  public :
  template<typename GfType> high_frequency_expansion_view(GfType const & x): base_type(x) {};
  template<typename F> void set_from_function(F f) { base_type::set_from_function(f);} // bug of autodetection in triqs::lazy on gcc   
  template<typename RHS> high_frequency_expansion_view & operator = (RHS const & rhs) { base_type::operator = (rhs); return *this; } 
  std::ostream & print_for_lazy(std::ostream & out) const { return out<<"high_frequency_expansion_view";}
  //int order_min() const { return this->domain().order_min();} 
  //int order_max() const { return this->domain().order_max();} 
 };

 class high_frequency_expansion : public impl::gf_impl<meshes::tail,false> {
  typedef impl::gf_impl<meshes::tail,false> base_type;
  public : 
  high_frequency_expansion():base_type() {} 
  high_frequency_expansion (size_t N1, size_t N2, base_type::domain_type const & domain_, base_type::indices_type const & indices_) : 
   base_type(N1,N2,domain_,indices_) {}
  template<typename GfType> high_frequency_expansion(GfType const & x): base_type(x) {}
  template<typename RHS> high_frequency_expansion & operator = (RHS const & rhs) { base_type::operator = (rhs); return *this; }
  //int order_min() const { return this->domain().order_min();} 
  //int order_max() const { return this->domain().order_max();} 
 };
*/

 // a trait to identity the tail
 template <typename G> struct is_tail : boost::is_base_of< tag::is_gf<meshes::tail>, G> {};

 typedef tqa::matrix_view<std::complex<double> , tqa::Option::Fortran> tail_result_type;

 /* -------------------------------------------
  *  Structure of algebra for algebra valued functions
  * ------------------------------------------ */
 struct tail_algebra_function_desc { 

  template<typename S> struct scalar {
   S s; scalar( S const & x) : s(x) {}
   template <typename T> struct call_rtype { typedef S type; };
   template<typename T> typename call_rtype<T>::type operator() (T) const {return s;}
   meshes::tail domain () const { return meshes::tail();}
  };

  template<typename ProtoTag, typename L, typename R> struct binary_node  { 
   L const & l; R const & r; binary_node (L const & l_, R const & r_):l(l_),r(r_) {}
   template <typename T> struct call_rtype {
    typedef tup::_binary_ops_<ProtoTag, typename tup::call_result_type<L,T>::type , typename tup::call_result_type<R,T>::type  > ops_type;
    typedef typename ops_type::result_type type;
   };
   template<typename T> typename call_rtype<T>::type operator() (T const & arg) const {return call_rtype<T>::ops_type::invoke(l(arg),r(arg));}
   meshes::tail domain () const { return l.domain();} // assert here domain are equal
  }; 

  template<typename L> struct negate  { 
   L const & l; negate (L const & l_):l(l_) {} 
   template <typename T> struct call_rtype {
    typedef tup::_unary_ops_<p_tag::negate, typename tup::call_result_type<L,T>::type > ops_type;
    typedef typename ops_type::result_type type;
   };
   template<typename T> typename call_rtype<T>::type operator() (T const & arg) const {return call_rtype<T>::ops_type::invoke(l(arg));}
   meshes::tail domain () const { return l.domain();}
  };
  
 };

 // specialize the * operation....
 template< typename L,typename R> struct tail_algebra_function_desc::binary_node<p_tag::multiplies,L,R>  { 
  L const & l; R const & r; binary_node (L const & l_, R const & r_):l(l_),r(r_) {}
  template <typename T> struct call_rtype { typedef tail_result_type type;};
  meshes::tail domain () const { 
   int omin = l.order_min()+r.order_min(); 
   return meshes::tail(omin,omin + std::min(l.domain().len(),r.domain().len()));
  } 
  tail_result_type operator() (int n) const {
   //if (N1!=t.N1 || N2!=t.N2) TRIQS_RUNTIME_ERROR<<"Multiplication is valid only for similar tail shapes !";
   //if (new_ordermin < OrderMinMIN) TRIQS_RUNTIME_ERROR<<"The multiplication makes the new tail have a too small OrderMin";
   tail_result_type::non_view_type res(l(l.domain().order_min()).shape()[0], l(l.domain().order_min()).shape()[1]); res()=0;
   const int kmin = std::max(0, n - r.domain().order_max() - l.domain().order_min() );
   const int kmax = std::min(l.domain().order_max() - l.domain().order_min(), n - r.domain().order_min() - l.domain().order_min() );
   std::cout<< kmin << "  "<<kmax<<std::endl;
   for (int k = kmin; k <= kmax; ++k)  { std::cout <<" k = "<< res<<std::endl; res += l(l.domain().order_min() +k) * r( n- l.domain().order_min() -k);}
   return res;
  }
 };

 template< typename L,typename R> struct tail_algebra_function_desc::binary_node<p_tag::divides,L,R>  { 
  L const & l; R const & r; binary_node (L const & l_, R const & r_):l(l_),r(r_) {}
  template<typename T> tail_result_type operator() (T const & arg) const {return l(arg)/r(arg);}
 };

 //TRIQS_PROTO_DEFINE_ALGEBRA_VALUED_FNT_ALG_WITH_DESC (is_tail, is_scalar_or_element, tail_algebra_function_desc);

 namespace impl { 
  template <typename Expr> struct tail_Expr;

  typedef triqs::utility::proto::algebra::grammar_generator<tail_algebra_function_desc,is_tail>::type grammar;
  typedef triqs::utility::proto::domain<grammar,tail_Expr,false>  tail_domain;

  template<typename Expr> struct tail_Expr : boost::proto::extends<Expr, tail_Expr<Expr>, tail_domain>{
   tail_Expr( Expr const & expr = Expr() ) : boost::proto::extends<Expr, tail_Expr<Expr>, tail_domain>  ( expr ) {}
   typedef typename boost::result_of<grammar(Expr) >::type _G;
   template<typename T> typename tup::call_result_type<_G,T>::type operator() (T const & x) const { return grammar()(*this)(x); }
   meshes::tail domain() const { return grammar()(*this).domain(); }
   friend std::ostream &operator <<(std::ostream &sout, tail_Expr<Expr> const &expr) { return boost::proto::eval(expr, triqs::utility::proto::AlgebraPrintCtx (sout)); }
   //int order_min() const { return this->domain().order_min();} 
   //int order_max() const { return this->domain().order_max();} 
  };
 }
 BOOST_PROTO_DEFINE_OPERATORS(is_tail, impl::tail_domain);


}}}

#endif

