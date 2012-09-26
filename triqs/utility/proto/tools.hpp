/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by O. Parcollet
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
#ifndef TRIQS_UTILITY_PROTO_H
#define TRIQS_UTILITY_PROTO_H

#include <boost/type_traits/add_const.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/proto/core.hpp>
#include <boost/proto/context.hpp>
#include <boost/proto/transform/arg.hpp>
#include <boost/proto/transform.hpp>
#include <boost/proto/domain.hpp>
#include <complex>
#include "triqs/utility/typeid_name.hpp"
#include <triqs/utility/view_tools.hpp>
#include <assert.h>

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>

namespace triqs { namespace utility { namespace proto { 
 namespace mpl = boost::mpl; namespace proto = boost::proto; namespace p_tag= proto::tag;
 namespace bf = boost::fusion;

 template<typename T> struct is_in_ZRC : boost::is_arithmetic<T>  {};
 template<> struct is_in_ZRC<bool> : mpl::true_ {};
 template<typename T> struct is_in_ZRC<std::complex<T> > :  mpl::true_ {};

 template <typename T> std::ostream & formal_print(std::ostream & out, T const & x) { return out<<x;}
 
 // to be moved in the genral lib  
 template<typename T> struct remove_const_and_ref : boost::remove_const<typename boost::remove_reference<T>::type> {};

 /* -------------------------------------------
  *   Print context
  * ------------------------------------------ */
 struct algebra_print_ctx : proto::callable_context< algebra_print_ctx const > {
  typedef std::ostream &result_type;
  result_type out;
  algebra_print_ctx(std::ostream & out_):out(out_) {}
  template <typename T>
   typename boost::enable_if<is_in_ZRC<T>, result_type>::type 
   operator ()(proto::tag::terminal, const T & A) const { return out<<A; }
  template <typename T>
   typename boost::disable_if<is_in_ZRC<T>, result_type>::type 
   operator ()(proto::tag::terminal, const T & A) const {return formal_print(out,A); }
  template<typename L, typename R> result_type operator ()(proto::tag::plus, L const &l, R const &r) const { return out << '(' << l << " + " << r << ')'; }
  template<typename L, typename R> result_type operator ()(proto::tag::minus, L const &l, R const &r) const { return out << '(' << l << " - " << r << ')'; }
  template<typename L, typename R> result_type operator ()(proto::tag::multiplies, L const &l, R const &r) const { return out << l << " * " << r; }
  template<typename L, typename R> result_type operator ()(proto::tag::divides, L const &l, R const &r) const { return out << l << " / " << r; }
 };

 template< typename T> std::ostream & print_algebra (std::ostream &sout, T const &x) {return proto::eval(x, algebra_print_ctx (sout));}

 /* ---------------------------------------------------------------------------------------------------
  * The domain can enforce copies or not...
  * --------------------------------------------------------------------------------------------------- */

 template< typename Grammar, template<typename Expr> class The_Expr> struct domain_with_ref: proto::domain<proto::generator<The_Expr>, Grammar> {};

 //If true it modifies the PROTO Domain to make *COPIES* of ALL objects.
 //cf http://www.boost.org/doc/libs/1_49_0/doc/html/proto/users_guide.html#boost_proto.users_guide.front_end.customizing_expressions_in_your_domain.per_domain_as_child
 //Objects which have a view type are however NOT copied, the copy is replaced by a VIEW.
 template< typename Grammar, template<typename Expr> class The_Expr> struct domain_with_copy : proto::domain<proto::generator<The_Expr>, Grammar> {
  //template< typename T > struct as_child : proto::domain<proto::generator<The_Expr>, Grammar>::proto_base_domain::template as_expr< T > {};
  template< typename T > struct as_child : 
   domain_with_copy<Grammar,The_Expr>::proto_base_domain::template as_expr< typename boost::add_const<typename const_view_type_if_exists_else_type <T>::type >::type > {};
 };

 // compat . to be removed 
 template< typename Grammar, template<typename Expr> class The_Expr, bool CopyOrViewInsteadOfRef> struct domain;

 template< typename Grammar, template<typename Expr> class The_Expr> struct domain<Grammar,The_Expr,false> : proto::domain<proto::generator<The_Expr>, Grammar> { };
 //If true it modifies the PROTO Domain to make *COPIES* of ALL objects.
 //cf http://www.boost.org/doc/libs/1_49_0/doc/html/proto/users_guide.html#boost_proto.users_guide.front_end.customizing_expressions_in_your_domain.per_domain_as_child
 //Objects which have a view type are however NOT copied, the copy is replaced by a VIEW.
 template< typename Grammar, template<typename Expr> class The_Expr> struct domain<Grammar,The_Expr,true> : proto::domain<proto::generator<The_Expr>, Grammar> {
  //template< typename T > struct as_child : proto::domain<proto::generator<The_Expr>, Grammar>::proto_base_domain::template as_expr< T > {};
  template< typename T > struct as_child : 
   domain<Grammar,The_Expr,true>::proto_base_domain::template as_expr< typename boost::add_const<typename const_view_type_if_exists_else_type <T>::type >::type > {};
 };


 /* ---------------------------------------------------------------------------------------------------
  * Simple transform to make * and /
  * --------------------------------------------------------------------------------------------------- */

 // a trick to avoid putting T() in the typeof type deduction ! This code is NEVER used. It is a hack to avoid decltype....
 template<class T> typename boost::remove_reference<typename boost::unwrap_reference<T>::type>::type * pseudo_default_construct() { 
  typename boost::remove_reference<typename boost::unwrap_reference<T>::type>::type * x= NULL; assert(0); return x; 
 }

 template<typename A,  typename B> 
  struct type_of_mult{
   typedef typename boost::remove_reference<A>::type T1; typedef typename boost::remove_reference<B>::type T2;
   typedef BOOST_TYPEOF_TPL( (*(pseudo_default_construct<T1>())) * (*(pseudo_default_construct<T2>()))) type;
  };

 struct negate_t { 
  BOOST_PROTO_CALLABLE();
  template<typename Sig> struct result;
  template<typename This, typename A> struct result<This(A)> : boost::remove_reference<A> {};
  template<typename A> typename result<negate_t(A)>::type  operator ()(A const & a) const { return - a; }
 };
 
 struct multiplies_t { 
  BOOST_PROTO_CALLABLE();
  template<typename Sig> struct result;
  template<typename This, typename A,  typename B> struct result<This(A,B)> : type_of_mult<A,B> {}; 
  template<typename A, typename B> typename result<multiplies_t(A,B)>::type  operator ()(A const & a, B const & b) const { return a*b; }
 };

 struct divides_t { 
  BOOST_PROTO_CALLABLE();
  template<typename Sig> struct result;
  template<typename This, typename A,  typename B> struct result<This(A,B)> { 
   typedef typename boost::remove_reference<A>::type T1; typedef typename boost::remove_reference<B>::type T2;
   typedef BOOST_TYPEOF_TPL( T1() * T2()) type;
  };
  template<typename A, typename B> typename result<divides_t(A,B)>::type  operator ()(A const & a, B const & b) const { return a/b; }
 };

 /* ---------------------------------------------------------------------------------------------------
  * A transform to evaluate ...
  * --------------------------------------------------------------------------------------------------- */

 /*
 // transform to eval the function call 
#define AUX6(z,p,unused) (*(pseudo_default_construct<typename boost::remove_reference<A##p>::type>()))
#define AUX(z, NN, unused) \
struct call_fnt_t##NN {\
BOOST_PROTO_CALLABLE();\
template<typename Sig> struct result;\
template<typename This, typename F BOOST_PP_ENUM_TRAILING_PARAMS(NN, typename A)> struct result<This(F BOOST_PP_ENUM_TRAILING_PARAMS(NN,A))> {\
typedef BOOST_TYPEOF_TPL ( (*(pseudo_default_construct<const F>())) (BOOST_PP_ENUM(NN,AUX6,nil))) type;\
};\
template<typename F BOOST_PP_ENUM_TRAILING_PARAMS(NN, typename A)>\
typename result<call_fnt_t##NN(F BOOST_PP_ENUM_TRAILING_PARAMS(NN,A))>::type\
operator()(const F & f BOOST_PP_ENUM_TRAILING_BINARY_PARAMS(NN,A, const &a)) const { return f(BOOST_PP_ENUM_PARAMS(NN,a));}\
};
BOOST_PP_REPEAT(BOOST_PP_INC(TRIQS_LAZY_MAXNARGS_CALLABLE), AUX, nil)
#undef AUX
#undef AUX6
*/

template<int Arity> struct eval_fnt;

#define AUX6(z,p,unused) (*(pseudo_default_construct<typename boost::remove_reference<A##p>::type>()))
template<> struct eval_fnt<1> {
 BOOST_PROTO_CALLABLE();
 template<typename Sig> struct result;
 //template<typename This, typename F, typename AL> struct result<This(F,AL)>
 // : boost::result_of<typename boost::remove_reference<F>::type 
 // (typename boost::remove_reference<typename bf::result_of::at< typename boost::remove_reference<AL>::type , mpl::int_<0> >::type>::type )> {};
 template<typename This, typename F, typename AL> struct result<This(F,AL)> {
  typedef typename bf::result_of::at< typename boost::remove_reference<AL>::type , mpl::int_<0> >::type A0;
  typedef BOOST_TYPEOF_TPL ( (*(pseudo_default_construct<const F>())) (BOOST_PP_ENUM(1,AUX6,nil))) type;
 };

 template<typename F, typename AL>
  typename result<eval_fnt(F,AL)>::type 
  operator ()(F const &f, AL const & al) const { return f(bf::at<mpl::int_<0> >(al)); }
};
#undef AUX6

 // Debug tool
 template<typename Expr>
  std::ostream & print_structure (std::ostream &sout, Expr const & E)  { 
   sout<<"Expression  "<<E <<std::endl;//triqs::utility::typeid_name(*this)<<std::endl;
   sout<<"            : return_type : "<<triqs::utility::typeid_name((typename Expr::value_type*)0)<<std::endl;
   sout<<"            : domain_type : "<<triqs::utility::typeid_name((typename Expr::domain_type*)0)<<std::endl;
   return sout;
  }

}}}

#endif

