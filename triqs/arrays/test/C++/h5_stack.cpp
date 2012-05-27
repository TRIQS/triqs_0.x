
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
#include "./src/h5/array_stack.hpp"
#include "./src/h5/simple_read_write.hpp"
#include <iostream>
#include "./src/impl/asserts.hpp"

using namespace std;
using namespace triqs::arrays;

template < class T>
void test(std::string filename, T init) { 

 h5::H5File file( filename.c_str(), H5F_ACC_TRUNC );

 const size_t N = 12, bufsize = 5, d= 2;

 array<T,2> A(d,d+1);
 array<T,3> A_stack_keep(N,d,d+1), A_stack_compare(N,d,d+1);

 array<T,1> B(d);
 array<T,2> B_stack_keep(N,d), B_stack_compare(N,d);

 T C;
 array<T,1> C_stack_keep(N), C_stack_compare(N);

 //h5::array_stack<T,2> SA( file, "A", A.shape() , bufsize);
 //h5::array_stack<T,1> SB( file, "B", B.shape() , bufsize);
 //h5::array_stack<T,0> SC( file, "C", mini_vector<size_t,0>() , bufsize);
 h5::array_stack< array<T,2 > > SA( file, "A", A.shape() , bufsize);
 h5::array_stack< array<T,1 > > SB( file, "B", B.shape() , bufsize);
 h5::array_stack< T> SC( file, "C", bufsize);
 //h5::array_stack< T> SC( file, "C", mini_vector<size_t,0>() , bufsize); // also valid ...
 //h5::array_stack< array<T,1 > > SB2( file, "B", bufsize); // does not compile
 
 for (int u = 0; u<N; ++u)  {
  A() = double(u+1)* init; 
  B() = double(u+1)* init; 
  C = double(u+1)* init; 
  A(0,0) *=2;
  B(0) *=2;
  SA() =  A; ++SA; 
  //SB() =  B; ++SB;
  SB << B;
  SC << C;
  A_stack_keep(u,range(),range()) = A;
  B_stack_keep(u,range()) = B;
  C_stack_keep(u) = C;
 }
 SA.flush();
 SB.flush();
 SC.flush();
 file.close(); // end writing 

 // now we read the file and compare

 h5::H5File file2( filename.c_str() ,H5F_ACC_RDONLY );
 h5_read (file2, "A",A_stack_compare); 
 cerr<< "A keep      "<< A_stack_keep<<endl;
 cerr<< "A in file   "<< A_stack_compare<<endl;
 
 assert_all_close (A_stack_keep, A_stack_compare, 1.e-13);

 h5_read (file2, "B",B_stack_compare); 
 cerr<< "B keep      "<< B_stack_keep<<endl;
 cerr<< "B in file   "<< B_stack_compare<<endl;
 
 assert_all_close (B_stack_keep, B_stack_compare, 1.e-13);

 h5_read (file2, "C",C_stack_compare); 
 cerr<< "C keep      "<< C_stack_keep<<endl;
 cerr<< "C in file   "<< C_stack_compare<<endl;
 
 assert_all_close (C_stack_keep, C_stack_compare, 1.e-13);

 file2.close();

}

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);

 test("stack_d.h5", 1.0 );
 test("stack_c.h5", complex<double>(1,2) );

}

