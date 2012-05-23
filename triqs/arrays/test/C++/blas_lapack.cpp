
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
#include "./src/linalg/matmul.hpp"
#include "./src/linalg/mat_vec_mul.hpp"
#include "./src/linalg/inverse.hpp"
#include "./src/linalg/determinant.hpp"

#include <boost/numeric/bindings/std/vector.hpp>
#include <boost/numeric/bindings/blas/level1/axpy.hpp>
#include <boost/numeric/bindings/blas/level2/ger.hpp>
#include <boost/numeric/bindings/blas/level2/gemv.hpp>
#include <boost/numeric/bindings/blas/level3/gemm.hpp>
#include <boost/numeric/bindings/lapack/computational/getrf.hpp>
#include <boost/numeric/bindings/lapack/computational/getri.hpp>


using namespace std;
using namespace triqs::arrays;

namespace blas = boost::numeric::bindings::blas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings= boost::numeric::bindings;

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);
 
 	typedef std::complex<double> T;
 triqs::arrays::vector<T>  V(5),V2(5);

 for (int i =0; i<5; ++i) {V(i) = i; V2(i) = -1;}

 cout<<"V = "<<V<<endl;
 cout<<"V2 = "<<V2<<endl;

 cout <<"starting blas test"<<endl;

 blas::axpy(2.0,V,V2);

 cout<<"V = "<<V<<endl;
 cout<<"V2 = "<<V2<<endl;

 triqs::arrays::vector <double> V3(2);
 for (int i =0; i<2; ++i) {V3(i) = i+1;}
 
 triqs::arrays::matrix<double,Option::Fortran > M1(2,2), M2(2,2), M3(2,2);
 for (int i =0; i<2; ++i)
  for (int j=0; j<2; ++j)
  { M1(i,j) = i+j; M2(i,j) = 1; M3(i,j)=0;}

 // try to multiply
 blas::gemm(1.0,M1, M2, 1.0, M3);
 cout<<"M1 = "<<M1<<endl;
 cout<<"M2 = "<<M2<<endl;
 cout<<"M3 = "<<M3<<endl;

/* triqs::arrays::matrix_view<double,'C' > MC1(M1);
 cout<<"M1 = "<<M1<<endl;
 cout<<"MC1 = "<<MC1<<endl;

*/

 // blas::get(alpha, A,B, c, R) : R <= alpha *A*B + c*R
 cout<<"V3 = "<<V3<<endl;
 blas::ger(1.0,V3,V3,M2);
 cout<<"M2 = "<<M2<<endl;

 // try to invert
 triqs::arrays::vector <int> ipiv(2);
 lapack::getrf(M1, ipiv);
 lapack::getri(M1, ipiv);
 cout<<"inv M1 = "<<M1<<endl;


}



