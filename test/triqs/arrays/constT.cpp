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
#include "./src/array.hpp"
#include <iostream>
#include <type_traits>
int main(int argc, char **argv) {

 
 triqs::arrays::array<long,2 > A (2,3);

 for (int i =0; i<2; ++i)
  for (int j=0; j<3; ++j) 
  { A(i,j) = 10*i+ j; }

 std::cerr<<"A = "<<A<<std::endl;
 
 triqs::arrays::array<const long,2 > B = A;

 A(0,0) +=20;
 //B = A; // DOES NOT AND SHOULD NOT COMPILE
 std::cerr<<"A = "<<A<<std::endl;
 
}




