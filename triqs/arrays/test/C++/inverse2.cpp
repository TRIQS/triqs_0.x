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
#include "./src/proto/matrix_algebra.hpp"
#include "./src/linalg/det_and_inverse.hpp"
#include "./src/linalg/matmul.hpp"
#include <iostream>

using namespace std;
using namespace triqs::arrays;
namespace blas = boost::numeric::bindings::blas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings= boost::numeric::bindings;

//using linalg::inverse;
//using linalg::inverse_and_compute_det;
//using linalg::determinant;

template<typename Expr >
matrix_view <typename Expr::value_type>
eval_as_matrix( Expr const & e) { return matrix<typename Expr::value_type>(e);}

int main(int argc, char **argv) {
 
 init_python_stuff(argc,argv);

 try { 

 triqs::arrays::matrix<double,Option::Fortran > W(3,3),Wi(3,3),Wkeep(3,3),A;
 for (int i =0; i<3; ++i)
  for (int j=0; j<3; ++j)
   W(i,j) = (i>j ? i+2.5*j : i*0.8-j);

 Wkeep = W;

 cout<< inverse(W) << " = "<< eval_as_matrix (inverse(W))<<endl<<endl;
 cout<< inverse(W) << " = "<< triqs::arrays::matrix<double,Option::Fortran >  (inverse(W))<<endl<<endl;
 cout<< inverse(W) << " = "<< triqs::arrays::matrix<double >  (inverse(W))<<endl<<endl;
 
 cout<< " and det = "<< double(determinant(W))<<endl<<endl;
 Wi = inverse(W);
 cout<< " Wi= "<< Wi<<endl<<endl;

 triqs::arrays::matrix<double,Option::Fortran > should_be_one(W*Wi);
  for (int i =0; i<3; ++i)
  for (int j=0; j<3; ++j)
    assert ( (abs(should_be_one(i,j) - (i==j ? 1 : 0))) <1.e-10 );

 cerr<< "W* inverse(W)" << " = "<< triqs::arrays::matrix<double,Option::Fortran > (W*Wi)<<endl<<endl;
 W=  inverse(W);
 cout<< " invert of W= "<< W<<endl<<endl;
 
 A=  inverse(W);
 cerr<< " A = inv W= "<< A<<endl<<endl;
 for (int i =0; i<3; ++i)
  for (int j=0; j<3; ++j)
    assert ( (abs(A(i,j) - Wkeep(i,j))) <1.e-10 );

 double det = determinant(W);
 cout<<"determinant "<<determinant(W)<< " = "<< det<< endl<<endl;

 //
 matrix_view<double,Option::Fortran> V(W(range(0,3,2), range(0,3,2)));
 cout<<" view = "<< V<<endl<<endl;
 cout<< inverse(V) << " = "<< eval_as_matrix (inverse(V))<<endl<<endl;


 // testing against "manual" call of bindings
 Wi = W;
 triqs::arrays::vector <int> ipiv2(3);
 lapack::getrf(Wi, ipiv2);
 cout<<"getrf W = "<<Wi<<endl<<endl;
 lapack::getri(Wi, ipiv2);
 cerr<<"inverse W = "<<Wi<<endl<<endl; // avoid printing because of 1.e-18 error, not reproducible
 for (int i =0; i<3; ++i)
  for (int j=0; j<3; ++j)
    assert ( (abs(Wi(i,j) - Wkeep(i,j))) <1.e-10 );


 }
 catch (std::string ERR) { cout<<"ERROR : "<< ERR<<endl;}

 return 0;

}



