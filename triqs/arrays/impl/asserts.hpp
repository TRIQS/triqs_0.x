
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

#ifndef TRIQS_ARRAYS_ASSERTS_H
#define TRIQS_ARRAYS_ASSERTS_H
#include "../expressions/array_algebra.hpp"
#include "../expressions/map.hpp"
#include "../expressions/min_max.hpp"

namespace triqs { namespace arrays {

 template<class T> inline double assert_abs(T z) { return std::abs(z);}

 template<class ArrayType1, class ArrayType2 >
  void assert_all_close( ArrayType1 const & A, ArrayType2 const & B, double precision) {
   typedef typename ArrayType1::value_type F;
   typename triqs::arrays::result_of::map<double (*)(F)>::type  Abs = map( static_cast< double (*)(F)> (assert_abs) );
   if ( max_element (Abs(A-B)) > precision) TRIQS_RUNTIME_ERROR<<"assert_all_close error : "<<A<<"\n"<<B;
  }

}}
#endif

