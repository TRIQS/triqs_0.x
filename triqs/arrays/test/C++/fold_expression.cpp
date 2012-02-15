
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

#include "./src/expressions/fold.hpp"
#include "./src/expressions/map.hpp"
#include "./src/expressions/matrix_algebra.hpp"
//#include "./src/expressions/arithmetic.hpp"
#include <iostream>
#include <functional>

using namespace std;
using namespace triqs::arrays;

#include "./src/expressions/min_max.hpp"
#include "./src/expressions/sum_prod.hpp"

// template<class T> T mmax(T const & x, T const &  y) { return std::max(x,y);}
//auto max_element = fold ( mmax<double>); 
// non C++0x version : 
// result_of::fold<double (*)(const double&, const double&)>::type max_element = fold ( mmax<double>); 

//auto Abs = map( static_cast< double (*)(double)> (std::abs) );
triqs::arrays::result_of::map<double (*)(double)>::type  Abs = map( static_cast< double (*)(double)> (std::abs) );

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);

 triqs::arrays::matrix<double,Option::Fortran > A(3,3),B(3,3),C;
 triqs::arrays::matrix<double > Ac(3,3);

 for (int i =0; i<3; ++i)
  for (int j=0; j<3; ++j)
  { A(i,j) = i+2*j+1; B(i,j) = i-3*j;}

 C = A+ B ;
 cout<< " A " << A<<endl;
 cout<< " B " << B<<endl;
 cout<< " A+B " << C<<endl;

 cout<< " max A : "<<triqs::arrays::max_element(A)<<endl;
 cout<< " max B : "<<triqs::arrays::max_element(B)<<endl;
 cout<< " max abs(B) : "<<triqs::arrays::max_element(Abs(B))<<endl;
 cout<< " max A+B : "<<triqs::arrays::max_element(A+B)<<endl;

 cout <<" sum(A) "<< sum(A)<<endl;
 cout <<" prod(A) "<< prod(A)<<endl;
 return 0;
}