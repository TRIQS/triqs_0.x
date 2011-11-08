
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
#include "./src/expressions/matrix_algebra.hpp"
#include <boost/type_traits/is_convertible.hpp>
#include <boost/proto/debug.hpp>
#include "./python_stuff.hpp"
#include <iostream>

using namespace std;
using namespace triqs::arrays;
using namespace indexmaps;
using namespace storages;

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);

  matrix <long> A (2,2), B(2,2),C(2,2);
  for (int i =0; i<2; ++i)
   for (int j=0; j<2; ++j) 
   { A(i,j) = 10*i+ j;B(i,j) = i+ 2*j;}
  cout<<" A = "<<A<<endl;
  cout<<" B = "<<B<<endl;

  cout<< (A+2*B) <<endl;
  cout<<"-------"<<endl;

  cout<< (A+2*B).domain()<<endl;
  cout<<"----EVAL ---"<<endl;

  C= A + 2*B;
  cout<<" C = A + 2*B = "<<C<<endl;

  C= std::plus<matrix<long> >()(A,B);
  cout<<" C = A+B =  "<<C<<endl;

  // array division
  matrix<double> Af (2,2), Bf(2,2),Cf(2,2);
  Af = A; Bf = B; Bf(0,0) = 1; Cf()=0;
  cout<<" Af = "<<Af<<endl;
  cout<<" Bf = "<<Bf<<endl;
 
  cout<< (Af* Bf).domain()<<endl;

  cout<<" computing  Cf = Af * Bf "<<endl;
  Cf = Af * Bf;
  cout<<" Cf  = "<<Cf<<endl;

  cout<<" matrix( Af * Bf )"<<matrix<double>(Af*Bf)<<endl;
  
  cout<<" matrix( Af * (Bf + Cf) )"<<matrix<double>(Af*(Bf+ Cf))<<endl;


 return 0;
}
