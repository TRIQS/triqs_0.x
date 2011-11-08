
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

#include "./src/vector.hpp"
#include <iostream>

using namespace std;
using namespace triqs::arrays;

int main(int argc, char **argv) {

 init_python_stuff(argc,argv);

 triqs::arrays::vector<double> V(5), W;

 for (int i =0; i<5; ++i) 
 V(i) = i+1;

 cout<< V<<endl; 
 W = V;
 cout<< W <<endl; 
 cout<< V<<endl; 

 triqs::arrays::vector<double> A(3);
 A() = 10; 
 W(range(0,6,2)) = A;
 cout<< W <<endl; 

 W +=V;
 cout<< W <<endl; 

 W -=V;
 cout<< W <<endl; 

 V *=2;
 cout<< "W = "<<W <<endl; 
 cout<< "V = "<<V <<endl; 
 triqs::arrays::swap(W,V);
 cout<< "W = "<<W <<endl; 
 cout<< "V = "<<V <<endl; 

 return 0;
}
