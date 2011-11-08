
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

#ifndef TRIQS_ARRAYS_FUNCTION_OBJECT_H
#define TRIQS_ARRAYS_FUNCTION_OBJECT_H

namespace triqs { namespace arrays { namespace function_object { 

 /**
  * Dressing a simple function into a function object
  */
 template<typename T> struct from_regular_function;
 template<typename R, typename A1> struct from_regular_function< R (*)(A1) > { 
  typedef R result_type;
  typedef A1 arg1_type;
  typedef A1 argument_type;
  static const unsigned int arity=1;
  R(*f) (A1);
  from_regular_function ( R(*f_)(A1) ) : f(f_) {}
  //R operator() (A1 const & a) const { return f(a);}
  R operator() (A1 a) const { return f(a);}
 };

 template<typename R, typename A1, typename A2> struct from_regular_function< R (*)(A1,A2) > { 
  typedef R result_type;
  typedef A1 arg1_type;
  typedef A2 arg2_type;
  typedef A1 first_argument_type;
  typedef A2 second_argument_type;
  static const unsigned int arity=2;
  R(*f) (A1,A2);
  from_regular_function ( R(*f_)(A1,A2) ) : f(f_) {}
  //R operator() (A1 const & a1, A2 const & a2) const { return f(a1,a2);}
  R operator() (A1 a1, A2 a2) const { return f(a1,a2);}
 };

 template<typename R, typename A1> 
  from_regular_function< R (*)(A1) > make_from_regular_function( R (*f)(A1)) { return  from_regular_function< R (*)(A1) >(f);} 
 template<typename R, typename A1, typename A2> 
  from_regular_function< R (*)(A1,A2) > make_from_regular_function( R (*f)(A1,A2)) { return  from_regular_function< R (*)(A1,A2) >(f);} 

 /**
  * Dressing a function template into a template struct that calls it and has the right result_type, etc......
  */
 template <typename ArgType, typename ReturnType, ReturnType(*F)(ArgType) > struct from_function_template {
  typedef ReturnType result_type;
  typedef ArgType arg_type;
  //result_type operator() (arg_type const & a) { return F(a);}
  result_type operator() (arg_type  a) { return F(a);}
 };

 /**
  * Dressing a class template into a template struct that calls it and has the right result_type, etc......
  */
 template <class T, template <class T2> class U > struct from_class_template { 
  typedef T result_type;
  typedef T arg_type;
  U<T> & f;
  from_class_template (U<T> const & f_): f(f_) {}
  //result_type operator() (arg_type const & a) { return f(a);}
  result_type operator() (arg_type a) { return f(a);}
 };

}}}//namespace triqs::arrays 

#endif

