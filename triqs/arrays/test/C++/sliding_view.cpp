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
#include "./src/sliding_view.hpp"
#include <iostream>
using namespace triqs::arrays;

int main(int argc, char **argv) {
 
 array<long,3> A (2,2,10), B;
 array<long,2> S(2,2); S()=0; S(0,0) = 123; S(1,1) = 45;

 A() = 0;
 for (int i =0; i<10; ++i) A(0,0,i) = i;
 
 B = A; A()=0;
 auto slv = make_sliding_view<2>(A);
 
 for (int i =0; i<10; ++i) { slv.set(i); slv(0,0) = i; }

 slv.set(0); 
 slv =125;
 slv = S;
 slv = S + 2* S;

 std::cerr  << " A = "<< A<< std::endl;
 std::cerr  << " B = "<< B<< std::endl;

 array_view<long,2> C(slv);
 std::cerr  << " C = "<< C.indexmap()<< std::endl;
 std::cerr  << " C = "<< C<< std::endl;
 }
