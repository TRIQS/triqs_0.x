
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
#include "./src/vector.hpp"
#include "./src/matrix.hpp"
#include "./src/linalg/eigenelements.hpp"

#include <iostream>
using namespace std;
using namespace triqs::arrays;
using namespace triqs::arrays::linalg;

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);

 matrix<double> A(3,3);

 for (int i =0; i<3; ++i)
  for (int j=0; j<=i; ++j)
  { A(i,j) = (i>j ? i+2*j : i-j); A(j,i) = A(i,j);}

 cerr<<"A = "<<A<<endl;
 eigenelements_worker< matrix_view <double>, true> w (A());
 w.invoke();
 cout<<"A = "<<A<<endl;
 cout<<" vectors = "<< w.values()<<endl;
 cout<<" values = "<< w.vectors()<<endl;

 A() =0;
 A(0,1) = 1;
 A(1,0) = 1;
 A(2,2) = 8;
 A(0,2) = 2;
 A(2,0) = 2;

 cout<<"A = "<<A<<endl;
 cout<<" values = "<< eigenelements(A(),true).first<<endl;
 cout<<" vectors = "<< eigenelements(A(),true).second<<endl;

 A() =0;
 A(0,1) = 1;
 A(1,0) = 1;
 A(2,2) = 8;

 cout<<"A = "<<A<<endl;
 cout<<" vectors = "<< eigenelements(A(),true).second<<endl;
 cout<<" values = "<< eigenelements(A(),true).first<<endl;
 cout<<"A = "<<A<<endl;

 return 0;

}



