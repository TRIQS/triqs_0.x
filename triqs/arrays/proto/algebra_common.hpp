
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_EXPRESSION_ALGEBRA_COMMON_H
#define TRIQS_ARRAYS_EXPRESSION_ALGEBRA_COMMON_H

#include <boost/mpl/int.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/proto/core.hpp>
#include <boost/proto/context.hpp>
#include <boost/proto/transform/arg.hpp>
#include <boost/proto/transform.hpp>
#include "../../utility/typeid_name.hpp"
#include <complex> 
#include "../impl/common.hpp"
#include "../impl/mini_vector.hpp"

namespace triqs { namespace arrays { namespace expressions { 

 namespace mpl = boost::mpl; namespace proto = boost::proto;
 using boost::enable_if; using  boost::enable_if_c;
 using proto::_left; using proto::_right; namespace p_tag= proto::tag;

 struct dC_ScalarGrammar : 
  proto::or_<
  proto::terminal< int >
  , proto::terminal< long >
  , proto::terminal< float >
  , proto::terminal< double >
  , proto::terminal< std::complex<float> >
  , proto::terminal< std::complex<double> >
  >
 {};

 /* -------------------------------------------
  *   Print context
  * ------------------------------------------ */
 template <typename T> std::ostream & formal_print(std::ostream & out, T const & x) { return out<<x;}

 struct AlgebraPrintCtx : proto::callable_context< AlgebraPrintCtx const > {
  typedef std::ostream &result_type;
  result_type out;
  AlgebraPrintCtx(std::ostream & out_):out(out_) {}
  template <typename T> result_type operator ()(proto::tag::terminal, const T & A) const { return formal_print(out,A); }
  template<typename L, typename R>
   result_type operator ()(proto::tag::plus, L const &l, R const &r) const { return out << '(' << l << " + " << r << ')'; }
  template<typename L, typename R>
   result_type operator ()(proto::tag::minus, L const &l, R const &r) const { return out << '(' << l << " - " << r << ')'; }
  template<typename L, typename R>
   result_type operator ()(proto::tag::multiplies, L const &l, R const &r) const { return out << l << " * " << r; }
  template<typename L, typename R>
   result_type operator ()(proto::tag::divides, L const &l, R const &r) const { return out << l << " / " << r; }
 };

  // Debug tool
 template<typename Expr>
  std::ostream & print_structure (std::ostream &sout, Expr const & E)  { 
   sout<<"Expression  "<<E <<std::endl;//triqs::utility::typeid_name(*this)<<std::endl;
   sout<<"            : return_type : "<<triqs::utility::typeid_name((typename Expr::value_type*)0)<<std::endl;
   sout<<"            : domain_type : "<<triqs::utility::typeid_name((typename Expr::domain_type*)0)<<std::endl;
   return sout;
  }
}}}//namespace triqs::arrays 
#endif

