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
#include "../proto/array_algebra.hpp"
#include "../functional/map.hpp"
#include "../algorithms.hpp"

namespace triqs { namespace arrays {

 template<class T> inline double assert_abs(T z) { return std::abs(z);}

 template<class ArrayType1, class ArrayType2 >
  void assert_all_close( ArrayType1 const & A, ArrayType2 const & B, double precision, bool relative = false) {
   typedef typename ArrayType1::value_type F;
   BOOST_AUTO(  Abs , map( boost::function<double(F)> (assert_abs<F>) ));
   auto r =  max_element (Abs(A-B));
   auto r2 =  max_element (Abs(A) + Abs(B));
   if ( r > (relative ? precision * r2 : precision) ) 
    TRIQS_RUNTIME_ERROR<<"assert_all_close error : \n\n"<<".. A = "<<A<<"\n\n"<<".. B= "<<B<<"\n\n"<< ".. Residue is r = "<<r;
  }

}}
#endif

