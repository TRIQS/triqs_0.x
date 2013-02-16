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
#ifndef TRIQS_CLEF_ARRAY_ADAPTER_H 
#define TRIQS_CLEF_ARRAY_ADAPTER_H
#include "triqs/clef2/clef.hpp"
#include "triqs/arrays/array.hpp"
#include "triqs/arrays/proto/array_algebra.hpp"

namespace triqs { 
 namespace clef { 

  namespace details { 
   template<int N> struct apply_function_on_minivector;
#define __ARG(z,N,unused) Arg[N]
#define IMPL(z, N, unused) \
   template<> struct apply_function_on_minivector<N> {\
    template<typename R, typename F, typename T>\
    static void invoke(R& r, F const &f,  arrays::mini_vector<T,N> const & Arg) { r = f(BOOST_PP_ENUM(N,__ARG,nil));}\
   };
   BOOST_PP_REPEAT_FROM_TO(1,BOOST_PP_INC(ARRAY_NRANK_MAX) , IMPL, nil);
#undef __ARG
#undef IMPL 
  }

  // redo using auto_assign !
  template <typename LHS, typename RHS>  // make arrays assignable 
   struct assignment_impl<true,LHS, RHS, typename boost::enable_if< triqs::arrays::is_amv_value_or_view_class<LHS> >::type > { 

    struct f { 
     RHS const & rhs;
     f(RHS const & rhs_):rhs(rhs_){}
     typedef typename LHS::value_type value_type;
     //template<typename KeyType> // convert to std::ptrdiff_t if needed ... normally not.  
     void operator()(value_type & p, arrays::mini_vector<std::ptrdiff_t,LHS::rank> const & key) const { 
      details::apply_function_on_minivector<LHS::rank>::invoke(p,rhs,key);}
    };

    static void invoke(LHS & lhs, RHS const & rhs) { triqs::arrays::indexmaps::foreach(f(rhs),lhs); };
   };
 } 
 namespace arrays { 

  template <typename T> 
   typename boost::enable_if<triqs::arrays::is_amv_value_or_view_class<T>, std::ostream>::type & 
   triqs_nvl_formal_print(std::ostream & out, T const & x) { return out<<x;}

 }
}

#endif

