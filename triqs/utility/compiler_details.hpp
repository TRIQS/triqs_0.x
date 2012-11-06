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
#ifndef TRIQS_COMPILER_DETAILS_H
#define TRIQS_COMPILER_DETAILS_H

#ifdef __GNUC__
#define restrict __restrict__ 
#endif

#ifdef __GXX_EXPERIMENTAL_CXX0X__
#define TRIQS_USE_STATIC_ASSERT
//#define USE_VARIADIC_TEMPLATES
#endif

#ifndef TRIQS_USE_STATIC_ASSERT
#include "boost/static_assert.hpp"
#define static_assert(X,MESS) BOOST_STATIC_ASSERT((X)) 
#endif

#endif

