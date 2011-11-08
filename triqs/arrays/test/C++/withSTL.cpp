
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
#include <iostream>

using namespace std;
using namespace triqs::arrays;

#include <vector>
#include <map>
#include <algorithm>

bool te(int x) { return (x<25);}

int main(int argc, char **argv) {
 init_python_stuff(argc,argv);

 array<long,2> A (2,3);
 for (int i =0; i<2; ++i)
  for (int j=0; j<3; ++j) 
   A(i,j) = 10*i+ j;

 std::vector<array<long,2> > VV;
 VV.push_back(A);

 map<string, array<long,2> > MAP;
 MAP["1"] = A;
 cout<<"printing map"<<endl;
 cout<<MAP["1"]<<endl;

 // Trying to put a vector in an array
 std::vector<int> V (10);
 array<int,1 > B(V.size()), C(V.size());
 
 //put something in V
 for (unsigned int i =0; i<10; ++i) V[i] = 10+i;

 // copy to B. Iterators on array are STL compliant so STL algorithms work.
 std::copy(V.begin(),V.end(),B.begin()); 
 cout<<" B = "<< B<<endl;

 // change B
 for (int i =0; i<10; ++i) B(i)*=2;
 
 // come back !
 std::copy(B.begin(),B.end(),V.begin()); 

 // print V
 std::copy(V.begin(),V.end(),std::ostream_iterator<int>(std::cout,","));

 // countif
 cout<<" B= "<< B<<endl;
 cout<<" Number of elements <25 : "<< std::count_if(B.begin(), B.end(),te)<<endl;
 //cout<<" Number of elements <25"<< std::count_if(B.begin(), B.end(),[](int x){ return x<25;});
 
 // max_element
 B(9) = 0; B(8) = 18;
 cout<<" B= "<< B<<endl;
 cout<<" max(B) "<< *std::max_element(B.begin(),B.end())<<endl;
 cout<<" max(B) position "<< std::max_element(B.begin(),B.end()).indices()[0]<<endl;

 // replace_if
 std::replace_if (B.begin(), B.end(), te, 0);
 cout<<" B= "<< B<<endl;
 
 //swap
 C()=0;
 cout<<" B= "<< B<<endl;
 cout<<" C= "<< C<<endl;
 std::swap(B,C);
 cout<<" B= "<< B<<endl;
 cout<<" C= "<< C<<endl;

 return 0;
}

