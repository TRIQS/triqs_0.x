
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
#include "./src/matrix.hpp"
#include "./src/linalg/matmul.hpp"
#include <iostream>

using std::cout; using std::endl;
using namespace triqs::arrays;

template<typename O1, typename O2, typename O3> void test(bool all =false)  { 
 matrix<double,O1> M1(2,2);
 matrix<double,O2> M2(2,2);
 matrix<double,O3> M3;
 for (int i =0; i<2; ++i)
  for (int j=0; j<2; ++j)
  { M1(i,j) = i+j; M2(i,j) = 1 + i -j ; }

 // The central instruction : note that matmul returns a lazy object 
 // that has ImmutableArray interface, and defines a specialized version assignment
 // As a result this is equivalent to matmul_with_lapack(M1,M2,M3) : there is NO intermediate copy.
 M3 = matmul(M1,M2);

 if (all) { 
  std::cout<<"M1 = "<<M1<<std::endl;
  std::cout<<"M2 = "<<M2<<std::endl;
  std::cout<<"M3 = "<<M3<<std::endl;
  std::cout<<"M4 = "<< matrix<double,Option::Fortran>(matmul(M1,M2)) <<std::endl;
  std::cout<<"M5 = "<< matrix<double>(matmul(M1,M2)) <<std::endl;


  for (int i =0; i<2; ++i)
   for (int j=0; j<2; ++j)
    M3(i,j) = matmul(M1,M2)[mini_vector<int,2>(i,j)];
 }

 std::cout<<"M3 = "<<M3<<std::endl<<"----------------"<<std::endl;

}

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);

 using namespace Option;
 test<C,C,C>(true);
 test<C,C,Fortran>();
 test<C,Fortran,Fortran>();
 test<C,Fortran,C>();
 test<Fortran,Fortran,Fortran>();
 test<Fortran,C,Fortran>();
 test<Fortran,Fortran,C>();
 test<Fortran,C,C>();
 return 0;
}
