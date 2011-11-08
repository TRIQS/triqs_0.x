
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

#ifndef TRIQS_ARRAYS_EXPRESSION_FOLD_H
#define TRIQS_ARRAYS_EXPRESSION_FOLD_H 
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include "./function_objects.hpp"
#include "../array.hpp"

namespace triqs { namespace arrays {

 template<class F> 
  class fold_worker  { 
   F f;
   typedef typename boost::remove_const<typename boost::remove_reference<typename F::result_type >::type>::type result_type;
   
   template<class A>
    struct fold_func_adaptor { 
     F f; result_type r;
     fold_func_adaptor(F f_, A a_):f(f_) { r= a_;}
     template<class KT> void operator()(A const & b, KT &) { r = f(r,b);}
    };

   public:

   fold_worker ( F const & f_):f(f_) {} 

   template<class A>   
    result_type operator() (A const & a, typename A::value_type init = typename A::value_type() )  const { 
     fold_func_adaptor<typename A::value_type> func(f,init);
     foreach(boost::ref(func),a);
     return func.r;
    }
  };

 template<class F> fold_worker<F> fold (F const & f) { return fold_worker<F>(f);}

 template<typename R, typename A1> 
  fold_worker< function_object::from_regular_function< R (*)(A1,A1) > >
  fold (R (*f)(A1,A1)) { return fold( function_object::make_from_regular_function(f)); }

 namespace result_of { 
  template<class F> struct fold { typedef fold_worker<F> type;};
 template<typename R, typename A1> 
   struct fold<R (*)(A1,A1) > { typedef  fold_worker< function_object::from_regular_function< R (*)(A1,A1) > > type;};
 }

}}//namespace triqs::arrays

#endif

