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
#include "./common.hpp"
#include <iostream>
#include <algorithm>

#include "./src/array.hpp"
#include "./src/vector.hpp"
#include "./src/matrix.hpp"
#include "./src/asserts.hpp"

using namespace triqs::arrays;

int main(int argc, char **argv) {

#ifndef TRIQS_WORKAROUND_INTEL_COMPILER_BUGS

 {
  vector<double> a = { 1,3,2}, b= {2,3,1};
  assert_is_true  ( std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end()));
  assert_is_false ( std::lexicographical_compare(b.begin(), b.end(), a.begin(), a.end()));
 }
 {
  vector<double> a = {1,3,2}, b= {1,2,10};
  assert_is_false  ( std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end()));
  assert_is_true   ( std::lexicographical_compare(b.begin(), b.end(), a.begin(), a.end()));
 }

#endif

}



