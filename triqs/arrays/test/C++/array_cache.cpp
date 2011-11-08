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

#include "./src/array.hpp"
#define TRIQS_ARRAYS_CACHE_VERBOSE
#include "./src/impl/cache.hpp"
#include "./src/expressions/array_algebra.hpp"
#include "./python_stuff.hpp"
#include <iostream>

using namespace std;
using namespace triqs::arrays;

void f( array_view<long,2> & A, long u) { 
 A(0,0) = u;
}

using details::cache_CC;
using details::cache_CR;
using details::cache_FC;
using details::cache_FR;

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);
 array<long,2, Option::Fortran > A (2,3);

 for (int i =0; i<2; ++i)
  for (int j=0; j<3; ++j) 
   A(i,j) = 10*i+ j;

 for (int i =0; i<2; ++i)
  for (int j=0; j<3; ++j) 
   cout<<"A "<<A(i,j)<<endl;

 array<long,2 > Ac (A);

 cout<<"A = "<<A<<endl;
 cout<<"Ac = "<<Ac<<endl;

 cout<<"----------"<<endl;
 cout<<"Expr = "<< cache_CC(2*A).view() <<endl;

 cout<<"----------"<<endl;
 cout<<"A F->C = "<< cache_CC(A).view() <<endl;

 cout<<"----------"<<endl;
 cout<<"A C->C = "<< cache_CC(Ac).view() <<endl;

 cout<<"----------"<<endl;
 f(cache_CR(Ac),287);
 cout<<"Ac = "<<Ac<<endl;

 cout<<"----------"<<endl;
 f(cache_CR(A),287);
 cout<<"A = "<<A<<endl;

 cout<<"----------"<<endl;
 cout<<"A = "<<A<<endl;
 cout<<"A(range(0,1),range(1,2) ) = "<<A( range(0,1),range(1,2) )<<endl;

 array_view <long,2,Option::Fortran> V(A(range(0,1), range(1,2)));
 array_view <long,2> Vc(Ac(range(0,1), range(1,2)));

 Vc(0,0) = 300;

 V = Vc;
 cout<<"A = "<<A<<endl;
 cout<<"Ac = "<<Ac<<endl;

 f(cache_CR(V),156);
 cout<<"A = "<<A<<endl;
 cout<<"Ac = "<<Ac<<endl;


 return 0;
}




