
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
#include <iostream>

#include "./src/array.hpp"
#include "./src/vector.hpp"
#include "./src/matrix.hpp"
#include "./src/proto/matrix_algebra.hpp"
#include "./src/linalg/matmul.hpp"
#include "./src/linalg/mat_vec_mul.hpp"
#include "./src/linalg/inverse.hpp"
#include "./src/linalg/determinant.hpp"
#include "./src/linalg/a_x_ty.hpp"
#include <boost/numeric/bindings/blas/level1/dot.hpp>


using namespace std;
using namespace triqs::arrays;
using linalg::a_x_ty;

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);

 triqs::arrays::matrix<double, Option::Fortran > A(5,5);
 triqs::arrays::matrix<double > Ac(5,5);
 typedef triqs::arrays::vector<double> vector_type;
 vector_type  MC(5), MB(5);
 A()= 0;

 for (int i =0; i<5; ++i)
  {  MC(i) = i;MB(i)=10*(i+1); }

 range R(1,3);
 cout<<" MC(R) = "<<MC(R)<< endl<<endl;
 cout<<" MB(R) = "<<MB(R)<< endl<<endl;

 A(R,R) += a_x_ty(1.0,MC(R),MB(R));
 cout<<" A(R,R) = "<<A(R,R)<< endl<<endl;

 A(R,R) += a_x_ty(1.0,MB(R),MC(R));

 cout<<" A(R,R) = "<<A(R,R)<< endl<<endl;

 A(R,R) = a_x_ty(1.0,MB(R),MC(R));

 cout<<" A(R,R) = "<<A(R,R)<< endl<<endl;
 
 cout<<" full A"<< A<<endl<<endl;


 cout<< " MB, MC, dot "<< MB << MC << boost::numeric::bindings::blas::dot(MB,MC)<<endl;
 cout<< " MC, MC, dot "<< MB << MC << boost::numeric::bindings::blas::dot(MC,MC)<<endl;

}



