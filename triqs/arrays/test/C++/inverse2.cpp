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

using std::cout; using std::endl;
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

 std::cout<< inverse(W) << " = "<< eval_as_matrix (inverse(W))<<std::endl<<std::endl;
 std::cout<< inverse(W) << " = "<< triqs::arrays::matrix<double,Option::Fortran >  (inverse(W))<<std::endl<<std::endl;
 std::cout<< inverse(W) << " = "<< triqs::arrays::matrix<double >  (inverse(W))<<std::endl<<std::endl;
 
 std::cout<< " and det = "<< double(determinant(W))<<std::endl<<std::endl;
 Wi = inverse(W);
 std::cout<< " Wi= "<< Wi<<std::endl<<std::endl;

 triqs::arrays::matrix<double,Option::Fortran > should_be_one(W*Wi);
  for (int i =0; i<3; ++i)
  for (int j=0; j<3; ++j)
    assert ( (abs(should_be_one(i,j) - (i==j ? 1 : 0))) <1.e-10 );

 std::cerr<< "W* inverse(W)" << " = "<< triqs::arrays::matrix<double,Option::Fortran > (W*Wi)<<std::endl<<std::endl;
 W=  inverse(W);
 std::cout<< " invert of W= "<< W<<std::endl<<std::endl;
 
 A=  inverse(W);
 std::cerr<< " A = inv W= "<< A<<std::endl<<std::endl;
 for (int i =0; i<3; ++i)
  for (int j=0; j<3; ++j)
    assert ( (abs(A(i,j) - Wkeep(i,j))) <1.e-10 );

 double det = determinant(W);
 std::cout<<"determinant "<<determinant(W)<< " = "<< det<< std::endl<<std::endl;

 //
 matrix_view<double,Option::Fortran> V(W(range(0,3,2), range(0,3,2)));
 std::cout<<" view = "<< V<<std::endl<<std::endl;
 std::cout<< inverse(V) << " = "<< eval_as_matrix (inverse(V))<<std::endl<<std::endl;


 // testing against "manual" call of bindings
 Wi = W;
 triqs::arrays::vector <int> ipiv2(3);
 lapack::getrf(Wi, ipiv2);
 std::cout<<"getrf W = "<<Wi<<std::endl<<std::endl;
 lapack::getri(Wi, ipiv2);
 std::cerr<<"inverse W = "<<Wi<<std::endl<<std::endl; // avoid printing because of 1.e-18 error, not reproducible
 for (int i =0; i<3; ++i)
  for (int j=0; j<3; ++j)
    assert ( (abs(Wi(i,j) - Wkeep(i,j))) <1.e-10 );


 }
 catch (std::string ERR) { std::cout<<"ERROR : "<< ERR<<std::endl;}

 return 0;

}



