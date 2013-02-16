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
#ifndef TRIQS_STD_FUNCTION_CLEF_H
#define TRIQS_STD_FUNCTION_CLEF_H
#include <boost/function.hpp>

namespace triqs { namespace clef { 

 template<typename R, typename A1> struct is_callable_leaf <  boost::function<R(A1)> > : mpl::true_ {};
 template<typename R, typename A1, typename A2> struct is_callable_leaf <  boost::function<R(A1,A2)> > : mpl::true_ {};
 template<typename R, typename A1, typename A2, typename A3> struct is_callable_leaf <  boost::function<R(A1,A2,A3)> > : mpl::true_ {};

}}
#endif

