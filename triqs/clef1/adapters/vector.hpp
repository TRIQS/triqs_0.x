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
#ifndef TRIQS_CLEF_VECTOR_ADAPTER_H 
#define TRIQS_CLEF_VECTOR_ADAPTER_H
#include "triqs/clef/core.hpp"
#include <vector>

namespace triqs { namespace clef { 

 template<typename T> struct is_assignable_leaf<std::vector<T> > : mpl::true_{}; 
 template<typename T, typename F> struct has_auto_assign_subscript< std::vector<T>, F >: mpl::true_{}; 

 template<typename T, typename Fnt>
 void triqs_nvl_auto_assign_subscript (std::vector<T> & v, Fnt f) {
  for (size_t i=0; i<v.size(); ++i) triqs::clef::assignment_subscript(boost::ref(v[i]) , f(i));}
  //for (size_t i=0; i<v.size(); ++i) v[i] = f(i);}  //assignment(boost::ref(lhs[i]) , rhs(i));}

 template <typename T> std::ostream & triqs_nvl_formal_print(std::ostream & out, std::vector<T> const & x) { return out<<"vector";}
 
}}

#endif

