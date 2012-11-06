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
#ifndef TRIQS_CLEF_FUNCTION_H
#define TRIQS_CLEF_FUNCTION_H
#include "./core.hpp"
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

namespace triqs { namespace clef { 

 template<typename F> struct function;
#define AUX1(z, p, unused) X##p()
#define AUX(z, NN, unused) \
 template<typename ReturnType, BOOST_PP_ENUM_PARAMS(NN, typename T)>\
 struct function<ReturnType(BOOST_PP_ENUM_PARAMS(NN, T))>  {\
  typedef boost::function<ReturnType(BOOST_PP_ENUM_PARAMS(NN, T))> std_function_type;\
  template<typename Expr, BOOST_PP_ENUM_PARAMS(NN, typename X)>  \
  function(Expr const & _e, BOOST_PP_ENUM_PARAMS(NN, X)) : _exp(new Expr(_e)),_fnt_ptr(new std_function_type(make_function(_e BOOST_PP_ENUM_TRAILING(NN, AUX1,nil)))) {}\
  function():_fnt_ptr(new std_function_type ()){}\
  ReturnType operator()(BOOST_PP_ENUM_BINARY_PARAMS(NN, T, const & t)) const { return (*_fnt_ptr)(BOOST_PP_ENUM_PARAMS(NN, t));}\
  TRIQS_CLEF_HAS_AUTO_ASSIGN();\
  template<typename F> \
  friend void triqs_nvl_auto_assign (function & x, F _f) {\
   * (x._fnt_ptr) =  std_function_type (_f);\
   x._exp = boost::shared_ptr <void> (new typename F::expression_type (_f.expr));\
  }\
  TRIQS_CLEF_ADD_LAZY_CALL_WITH_COPY(NN, function);\
  friend std::ostream & triqs_nvl_formal_print(std::ostream & out, function const & x) { return out<<"function";}\
  private : \
	    boost::shared_ptr <void>  _exp;\
	     boost::shared_ptr < std_function_type > _fnt_ptr;\
 }; 

 BOOST_PP_REPEAT_FROM_TO(1,BOOST_PP_INC(TRIQS_CLEF_MAXNARGS_CALLABLE) , AUX, nil)
#undef AUX
#undef AUX1
  /*
   *   template<typename F> \
   void set_from_function (F _f) {\
   * _fnt_ptr =  std_function_type (_f);\
   _exp = boost::shared_ptr <void> (new typename F::expression_type (_f.expr));\
   }\
   */

}} //  namespace triqs::clef

#endif



