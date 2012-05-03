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

#define BOOST_RESULT_OF_USE_DECLTYPE
#include <boost/utility/result_of.hpp>

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
// return out<< triqs::utility::typeid_name(A); }
   template<typename L, typename R> result_type operator ()(proto::tag::plus, L const &l, R const &r) const { return out << '(' << l << " + " << r << ')'; }
   template<typename L, typename R> result_type operator ()(proto::tag::minus, L const &l, R const &r) const { return out << '(' << l << " - " << r << ')'; }
   template<typename L, typename R> result_type operator ()(proto::tag::multiplies, L const &l, R const &r) const { return out << l << " * " << r; }
   template<typename L, typename R> result_type operator ()(proto::tag::divides, L const &l, R const &r) const { return out << l << " / " << r; }
  };

 /* ---------------------------------------------------------------------------------------------------
  * The domain can enforce copies or not...
  * --------------------------------------------------------------------------------------------------- */

 template<typename T, typename Void =void> struct const_view_type_if_exists_else_type {typedef T type;}; 
 template<typename T> struct const_view_type_if_exists_else_type<T, typename T::has_view_type_tag> {typedef const typename T::view_type type;}; 
 // template<typename T> struct const_view_type_if_exists_else_type<const T, typename T::has_view_type_tag> {typedef const typename T::view_type type;}; 

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

  struct multiplies_t { 
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename A,  typename B> struct result<This(A,B)> { 
    typedef typename boost::remove_reference<A>::type T1; typedef typename boost::remove_reference<B>::type T2;
    typedef BOOST_TYPEOF_TPL( T1() * T2()) type;
   };
   template<typename A, typename B> typename result<multiplies_t(A,B)>::type  operator ()(A const & a, B const & b) const { return a*b; }
  };

  struct divides_t { 
   BOOST_PROTO_CALLABLE();
   template<typename Sig> struct result;
   template<typename This, typename A,  typename B> struct result<This(A,B)> { 
    typedef typename boost::remove_reference<A>::type T1; typedef typename boost::remove_reference<B>::type T2;
    typedef BOOST_TYPEOF_TPL( T1() / T2()) type;
   };
   template<typename A, typename B> typename result<divides_t(A,B)>::type  operator ()(A const & a, B const & b) const { return a/b; }
  };

  /* ---------------------------------------------------------------------------------------------------
  * A transform to evaluate ...
  * --------------------------------------------------------------------------------------------------- */

 template<int Arity> struct eval_fnt;
 template<> struct eval_fnt<1> {
  BOOST_PROTO_CALLABLE();
  template<typename Sig> struct result;
  template<typename This, typename F, typename AL> struct result<This(F,AL)>
   : boost::result_of<typename boost::remove_reference<F>::type 
   (typename bf::result_of::at< typename boost::remove_reference<AL>::type , mpl::int_<0> >::type )> {};

  template<typename F, typename AL>
   typename boost::result_of<F(typename bf::result_of::at<AL const, mpl::int_<0> >::type)>::type 
   operator ()(F const &f, AL const & al) const { return f(bf::at<mpl::int_<0> >(al)); }
 };


}}}

#endif

