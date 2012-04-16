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

namespace triqs { namespace gf { namespace local {

 namespace tqa= triqs::arrays; namespace tql= triqs::lazy; namespace mpl= boost::mpl;
 using tqa::range;

 namespace tag { struct is_gf{};}

 template<typename DomainType> class gf;
 template<typename DomainType> class gf_view;

 namespace impl {

  template<typename LHS, typename RHS, typename Enable=void > struct assignment;

  struct nothing { // used in place of tail when there is none... an object that does .. nothing 
   typedef nothing view_type;
   template<typename A, typename B> view_type slice(A a, B b) const { return nothing();}
   nothing(...){}
   template<typename T1, typename T2, typename T3, typename T4> nothing(T1,T2,T3,T4){}
   template<typename T> void operator=(T) {}
   template<typename T> void set_from_function(T) {}
  };

  template<typename TailViewType>
   struct wrap_infty : tag::is_gf {
    TailViewType tail;
    wrap_infty(TailViewType const & t): tail(t) {}
    //operator TailViewType() const { return tail; } 
    operator domains::infty() const { return domains::infty();}
    typename TailViewType::const_result_type operator() ( typename TailViewType::arg0_type const & x) const { return tail(x);}
   };

  /*------------------------------------------------------------------------------
   * The implementation class for both the view and the regular gf class
   *-----------------------------------------------------------------------------*/
  template<typename DomainType, bool IsView> class gf_impl : tag::is_gf {

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

    typedef typename mpl::if_c<domain_type::has_tail, gf<meshes::tail>, nothing>::type tail_non_view_type; 
    typedef typename tail_non_view_type::view_type                                     tail_view_type; 
    typedef typename mpl::if_c<IsView, tail_view_type, tail_non_view_type>::type       tail_type;

    typedef typename domain_type::index_type arg0_type;
    typedef tqa::matrix_view<value_element_type, arrays::Option::Fortran>       result_type;
    typedef tqa::matrix_view<const value_element_type, arrays::Option::Fortran> const_result_type;

   protected:
    domain_type _domain;
    data_type data;
    tail_type tail;
    indices_type _indices;

   public:

    domain_type const & domain() const { return _domain;}

    data_view_type data_view()             { return data;}
    const data_view_type data_view() const { return data;}

    tail_view_type tail_view()             { return tail;}
    const tail_view_type tail_view() const { return tail;}

    indices_type const & indices() const { return _indices;}

   protected:
    gf_impl (domain_type const & domain_, data_view_type const & data_, tail_view_type const & tail_, indices_type const & indices_) : 
     data(data_), _domain(domain_), tail(tail_), _indices(indices_){}

    gf_impl(gf_impl const & x): _domain(x._domain), data(x.data), tail(x.tail), _indices(x.indices()){}

    template<typename GfType> gf_impl(GfType const & x): _domain(x.domain()), data(x.data_view()), tail(x.tail_view()), _indices(x.indices()){}

    gf_impl() {} 

    gf_impl (size_t N1, size_t N2, domain_type const & domain_, indices_type const & indices_) : 
     _domain(domain_), data(N1,N2,domain().len()), tail(N1,N2,domain().mesh_tail,indices_), _indices(indices_){}

   public:

    triqs::arrays::mini_vector<size_t,2> result_shape() const { return data.shape();}

    result_type       operator() ( arg0_type const & x)       { return data(range(),range(),x);}
    const_result_type operator() ( arg0_type const & x) const { return data(range(),range(),x);}

    tail_view_type       operator() ( domains::infty const & x)       { return tail;}
    const tail_view_type operator() ( domains::infty const & x) const { return tail;}

    TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(1,view_type);

    template<typename A, typename B> view_type slice(A a, B b) { 
     return view_type( gf_impl<domain_type,true>(_domain, data (range(a), range(b), range()), tail.slice(a,b), slice_indices(indices(),range(a),range(b)) )); 
    }
    template<typename A, typename B> const view_type slice(A a, B b) const { 
     return view_type( gf_impl<domain_type,true>(_domain, data (range(a), range(b), range()), tail.slice(a,b), slice_indices(indices(),range(a),range(b)) )); 
    }

    template<typename RHS> void operator = (RHS const & rhs) { _domain = rhs.domain(); data = rhs.data; tail= rhs.data; } 

   private:
    template<typename F> void compute_tail(gf_impl<meshes::tail, IsView> & t, F f) { 
     t.set_from_function( f( wrap_infty<tail_view_type>( tail_non_view_type(data.shape()[0], data.shape()[1],domain().mesh_tail ,indices()))));
    }
    template<typename F> void compute_tail(nothing & t, F f) {}

   public:
    // lazy_assignable // TO DO : the computation of the tail is not ready
    template<typename F> void set_from_function(F f) { 
     const size_t Nmax = this->data.shape()[2]; for (size_t u=0; u<Nmax; ++u) this->data(range(),range(),u) = f(1.0*u);
     compute_tail(tail,f); // the tail is only computed if it exists, or f maybe called for infty in undefined cases.
    }

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


 // -------------------------------   Expression template --------------------------------------------------

 // a trait to identity the local gf
 template <typename G> struct is_gf : boost::is_base_of< tag::is_gf, G> {} ;
 //template <> struct is_gf<  impl::gf_impl<triqs::gf::meshes::matsubara_freq, true>::wrap_infty >: mpl::true_{};
 //static_assert( boost::is_base_of<tag::is_gf, impl::gf_impl<triqs::gf::meshes::matsubara_freq, true>::wrap_infty >::value, "oops");

 // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
 template <typename T> struct is_scalar_or_element : mpl::or_< tqa::expressions::matrix_algebra::IsMatrix<T>, triqs::utility::proto::is_in_ZRC<T> > {};

 namespace tupa=triqs::utility::proto::algebra;

 template <typename Expr> struct The_Expr;  // the expression
 typedef tupa::grammar_generator<tupa::algebra_function_desc,is_gf, is_scalar_or_element>::type grammar; // the grammar
 typedef tupa::domain<grammar,The_Expr,true>  domain; // the domain 

 template<typename Expr> struct The_Expr : boost::proto::extends<Expr, The_Expr<Expr>, domain>{ // impl the expression
  typedef boost::proto::extends<Expr, The_Expr<Expr>, domain> basetype;
  The_Expr( Expr const & expr = Expr() ) : basetype ( expr ) {}
  typedef typename boost::result_of<grammar(Expr) >::type _G;

  template<typename T> typename triqs::utility::proto::call_result_type<_G,size_t>::type operator() (T x) const { return grammar()(*this)(x); }

  // formal print of the expression
  friend std::ostream &operator <<(std::ostream &sout, The_Expr<Expr> const &expr) { return boost::proto::eval(expr, triqs::utility::proto::AlgebraPrintCtx (sout)); }
 };

 BOOST_PROTO_DEFINE_OPERATORS(is_gf, domain);

}}}

#endif

