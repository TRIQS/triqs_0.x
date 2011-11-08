
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

#ifndef TRIQS_ARRAYS_NUPLETOOLS_H
#define TRIQS_ARRAYS_NUPLETOOLS_H
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/int.hpp>
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"

namespace triqs { namespace arrays { namespace TupleTools { 

 /// Meta-function to count the number of elements of integral type
 template<typename TUPLE>
  struct CountHowManyInt { 
   static const int value = CountHowManyInt<typename TUPLE::tail_type>::value + 
    boost::mpl::if_<boost::is_integral<typename TUPLE::head_type>,boost::mpl::int_<1>,boost::mpl::int_<0> >::type::value;
  };

 template<> struct CountHowManyInt<boost::tuples::null_type > { static const int value=0; };

 /// Meta-function to count the number of elements of a given type
 template<typename T, typename TUPLE>
  struct CountHowMany { 
   static const int value = CountHowMany<T,typename TUPLE::tail_type>::value + 
    boost::mpl::if_<boost::is_base_of<T,typename TUPLE::head_type>,boost::mpl::int_<1>,boost::mpl::int_<0> >::type::value;
  };

 template<typename T> struct CountHowMany<T,boost::tuples::null_type > { static const int value=0; };


}}}//namespace triqs::arrays 
#endif
