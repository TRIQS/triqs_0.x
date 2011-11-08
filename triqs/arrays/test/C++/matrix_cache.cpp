
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

#include "./src/linalg/matrix_cache.hpp"
#include "./src/expressions/matrix_algebra.hpp"

#include <iostream>

using namespace std;
using namespace triqs::arrays;

template < class T>
void writeln(T const & x) { std::cout<<x<<std::endl;}

template <class T> 
void print (T const & x)   {

 matrix_cache<T,  matrix_cache_policy::const_copy > C(x); 
 std::cout<<C() <<std::endl;
}


int main(int argc, char **argv) {
 
 init_python_stuff(argc,argv);

 
 matrix<double> A (2,2), B(2,2);
for (int i =0; i<2; ++i)
  for (int j=0; j<2; ++j) 
   A(i,j) = 10*i+ j;

 {
  matrix_cache_details::matrix_cache< matrix<double>, matrix_cache_policy::const_copy>  C(A);
  writeln (C());
 }

{
 matrix_cache_details::matrix_cache< matrix<double>, matrix_cache_policy::reflexive>  C(A);
  writeln (C());
  C()(0,0) = 100;
 }
  writeln (A);

{
  matrix_cache< matrix_view<double>, matrix_cache_policy::reflexive>  C(A(range(0,1),range()));
  writeln (C());
  C()(0,0) = 33;
 }
  writeln (A);
{
  matrix_cache< matrix_view<double>, matrix_cache_policy::reflexive>  C(A(range(),range(0,1)));
  writeln (C());
  C()(0,0) = 34;
 }
  writeln (A);

  print (2*A);


 return 0;

}



