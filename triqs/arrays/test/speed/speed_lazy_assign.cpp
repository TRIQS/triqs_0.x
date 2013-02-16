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

using namespace triqs::arrays;
using namespace triqs;
const int N1= 200, N2 = 300;
namespace tqa=triqs::arrays; namespace tql=triqs::clef;

struct  aux_t { 
  double operator() (long i, long j) { return i + 2.0*j;}
 }; 

struct plain {
 void operator()() { 
  triqs::arrays::array<double,2> A (N1,N2,FORTRAN_LAYOUT);
  for (int u =0; u<5000; ++u)
  {
   for (int j=0; j<A.len(1); ++j) 
    for (int i =0; i<A.len(0); ++i)
     A(i,j) = i+ 2*j;
  }
 }
};


struct lazy {
 void operator()() {
  tql::placeholder<0> i_;   tql::placeholder<1> j_;  
  //triqs::arrays::array<double,2> A (N1,N2);
  triqs::arrays::array<double,2,TRAVERSAL_ORDER_FORTRAN> A (N1,N2,FORTRAN_LAYOUT);
  auto f = make_function(  i_+ 2.0*j_, i_, j_);
  aux_t aux;
  for (int u =0; u<5000; ++u)
   //indexmaps::foreach_av(boost::ref(aux), A);
   //indexmaps::foreach_av(f, A);
   A(i_,j_) << i_+ 2*j_;
 }
};

#include "./speed_tester.hpp"
int main() {
 const int l = 100;
 speed_tester<plain> (l);
 speed_tester<lazy> (l);
 return 0;
}

