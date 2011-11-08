
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
#include "./src/linalg/matmul.hpp"
#include <iostream>

using namespace std;
using namespace triqs::arrays;
namespace bindings= boost::numeric::bindings;

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);

 // testing gemv
 triqs::arrays::matrix<double,Option::Fortran > A(5,5);

 for (int i =0; i<5; ++i)
  for (int j=0; j<5; ++j)
   A(i,j) = i+2*j+1; 

 range R(1,3); 

 cout<< "A = "<< A<< endl;
 cout<<A.transpose() << endl;
 cout<<A(R,R) << endl;
 cout<<A(R,R).transpose() << endl;

 triqs::arrays::matrix_view<double > Acw =  A.transpose();
 cout<<"bindings::stride_major(a) "<<  bindings::stride_major(Acw(R,R)) <<endl; 
 cout<<"bindings::stride_major(a) "<<  bindings::stride_major(A(R,R)) <<endl; 
  
}



