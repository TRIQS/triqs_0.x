
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

#include "./src/expressions/map.hpp"
#include "./src/expressions/matrix_algebra.hpp"
//#include "./src/expressions/arithmetic.hpp"
#include <iostream>
#include <algorithm>

using namespace std;
using namespace triqs::arrays;
using namespace triqs::arrays::function_object;

template<class T> T mmax(T const & x, T const &  y) { return std::max(x,y);}

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);

 triqs::arrays::matrix<double,Option::Fortran > A(3,3),B(3,3);
 A() = -2;

 for (int i =0; i<3; ++i)
  for (int j=0; j<3; ++j)
  { A(i,j) = i+2*j+1; B(i,j) = i-j;}

 //auto Abs = map( static_cast< double (*)(double)> (std::abs) );
 triqs::arrays::result_of::map<double (*)(double)>::type  Abs = map( static_cast< double (*)(double)> (std::abs) );

 //auto Max = map( mmax<double> );
 triqs::arrays::result_of::map<double (*)(const double&, const double&)>::type Max = map( mmax<double> );

 cout<< " A " << A<<endl<<endl;
 cout<< " B " << B<<endl<<endl;
 cout<<" abs(B+B) = "<<eval(Abs(B+B)) <<endl<<endl;
 cout<<" A+10*B = "<<eval(A+10*B) <<endl<<endl;
 cout<<" Abs(A+10*B) = "<<eval(Abs(A+10*B)) <<endl<<endl;
 cout<<" Max(A,10*B)"<< eval(Max(A,10*B))<<endl<<endl;
 return 0;
}
