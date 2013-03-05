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
#include "./src/matrix.hpp"
#include "./src/array.hpp"
#include "./src/sliding_view.hpp"
using namespace std;
namespace tqa = triqs::arrays;

const int N = 1000;
const int nl_interne = 10000;

//typedef double VALUE_TYPE;
typedef int VALUE_TYPE;
inline VALUE_TYPE fnt(size_t i) { return i;} //*(i+2.0)*(i-8);}
//inline VALUE_TYPE fnt(size_t i) { return i*(i+2.0)*(i-8);}

struct array_code {
 void operator()() {

  tqa::array<VALUE_TYPE,3> A (2,2,N);
  A() = 0;
  for (int u =0; u<nl_interne; ++u)
   for (int i =0; i<N-1; ++i) A(0,0,i) = fnt(i);

 }
};

struct with_sliding_view {
 void operator()() {

  tqa::array<VALUE_TYPE,3> A (2,2,N);
  A() = 0;

  auto slv = tqa::make_sliding_view<2>(A);

  for (int u =0; u<nl_interne; ++u)
   for (int i =0; i<N-1; ++i ) {slv.set(i); slv(0,0) = fnt(i);}
 }
};
struct with_slices {
 void operator()() {

  tqa::array<VALUE_TYPE,3> A (2,2,N);
  A() = 0;

  for (int u =0; u<nl_interne; ++u)
  {
   for (int i =0; i<N-1; ++i) {A(tqa::range(),tqa::range(),i)(0,0) = fnt(i);}
  }
 }
};


#include "./speed_tester.hpp"
int main() {
 speed_tester<array_code> (5000);
 speed_tester<with_sliding_view> (5000);
 speed_tester<with_slices> (5000);
 return 0;
}

