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
#include "./python_stuff.hpp"
#include "./src/array.hpp"
#include <iostream>
#include <triqs/arrays/indexmaps/permutation2.hpp>
#include <triqs/arrays/indexmaps/cuboid/index_order.hpp>

using namespace triqs::arrays;
using namespace triqs::arrays::permutations;

int main(int argc, char **argv) {

  init_python_stuff(argc,argv);

 constexpr auto p=  permutation2(0,2,1);
 constexpr auto p2= permutation2(2,1,0,3);
 constexpr auto pc= permutation2(1,2,3,0);

 constexpr auto S1 = indexmaps::index_order::sliced_memory_order<int,int,range,int>(p2);

 std::cout  << " sliced "<< S1 << std::endl ; 
}
