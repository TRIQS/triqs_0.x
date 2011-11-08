
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

int main(int argc, char **argv) {
 init_python_stuff(argc,argv);

 //boost::function_requires<  boost::Mutable_ForwardIteratorConcept<array<long,2>::iterator > >();
 // boost::function_requires<  boost::ForwardIteratorConcept<array<long,2>::const_iterator > >();

 array<long,2> A (2,3);

 // first print the index generator
 for (array<long,2>::indexmap_type::domain_type::generator it = A.indexmap().domain().begin(); it; ++it)
  cout<<"  "<<*it<<endl;

 cout<<endl;

 for (indexmaps::cuboid_index_generator< array<long,2>::indexmap_type::domain_type, Permutations::permutation<0,1> >
   it(A.indexmap().domain()); it; ++it)
  cout<<"  "<<*it<<endl;

 cout<<endl;

 for (indexmaps::cuboid_index_generator< array<long,2>::indexmap_type::domain_type, Permutations::permutation<1,0> >
   it(A.indexmap().domain()); it; ++it)
  cout<<"  "<<*it<<endl;

 cout<<endl;

 cout<<" C order : traversal"<<endl;

 for (array<long,2>::iterator it = A.begin(); it; ++it) { 
  *it =it.indices()[0] + 10 *it.indices()[1] ;
  cout<<" i,j = "<<it.indices()<<endl;
 }
 cout <<"A = i + 10*j"<<A<<endl;

 int u=0;
 for (array<long,2>::iterator it = A.begin(); it; ++it,++u) { 
  *it =u;
  cout<<" i,j = "<<it.indices()<<endl;
 }
 cout <<"A = order of visit "<<A<<endl;


 cout<<" F order : traversal"<<endl; 
 array<long,2,Option::Fortran> Af (2,3);

 for (array<long,2,Option::Fortran>::iterator it = Af.begin(); it; ++it) { 
  *it =it.indices()[0] + 10 *it.indices()[1] ;
  cout<<" i,j = "<<it.indices()<<endl;
 }
 cout <<"A = i + 10*j"<<Af<<endl;

 u=0;
 for (array<long,2,Option::Fortran>::iterator it = Af.begin(); it; ++it,++u) { 
  *it = -u;
  cout<<" i,j = "<<it.indices()<<endl;
 }
 cout <<"Af = order of visit "<<Af<<endl;

 return 0;
}

