
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

#include "./src/expressions/function_objects.hpp"

/**
 * Examples
 */
#include <iostream>
using namespace std;
using namespace triqs::arrays::function_object;

// simple function 
double f(int i){ return i + 0.5;}

// template function
template <class T> T g(T x) { return x+1;} 
template<typename T> struct dressed_g : from_function_template<T,T,g> {};

// template class
template <class T> struct S { T operator()(T x) { return x+1;} };
template <class T> struct dressed_S : from_class_template<T,S> { dressed_S():from_class_template<T,S>(*this) {} };

int main() { 

 cout<< make_from_regular_function(f)(0)<<endl;
 cout<< dressed_g<double>() (2) <<endl;
 cout<< dressed_S<double>()(2)<<endl;

}
