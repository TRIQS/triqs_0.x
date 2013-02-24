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
#ifndef TRIQS_H5_VECTOR_H
#define TRIQS_H5_VECTOR_H
#include "./group.hpp"
#include <vector>

namespace triqs { namespace h5 {

 // to be done : vector of something else....

 // special case of vector of string
 inline void h5_write (group f, std::string const & name, std::vector<std::string> const & V) {
  detail::write_1darray_vector_of_string_impl(f,name,V);
 }

 inline void h5_read (group f, std::string const & name, std::vector<std::string> & V) {
  detail::read_1darray_vector_of_string_impl(f,name,V);
 }

}} 
#endif


