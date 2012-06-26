
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

#ifndef TRIQS_STORAGES_COMMON_H
#define TRIQS_STORAGES_COMMON_H
#include <boost/type_traits/remove_const.hpp>

namespace triqs { namespace arrays { 
 namespace Tag {struct storage{}; }
 namespace storages  { 

  // is a raw copy possible ? only if the value_type is the same
  template<typename S1, typename S2>
   struct raw_copy_possible : boost::is_same< typename boost::remove_const<typename S1::value_type>::type,
   typename boost::remove_const<typename S2::value_type>::type > {};

 }
}}//namespace triqs::arrays 
#endif

