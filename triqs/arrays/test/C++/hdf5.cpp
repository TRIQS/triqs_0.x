
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
#include "./src/h5/simple_read_write.hpp"

using namespace std;
using namespace triqs::arrays;

template <typename T> 
ostream & operator << (ostream & out, std::vector<T> const & v) { 
for (size_t i =0; i<v.size(); ++i) out<< v[i];
return out;
}

int main(int argc, char **argv) {
 
  init_python_stuff(argc,argv);

 try { 

 array<long,2> A (2,3),B,vc;
 array<double,2> D (2,3), D2;

 array<long,2,Option::Fortran> Af,Bf,vf;

 array<complex<double>,1> C(5), C2;
 complex<double> z(1,2);

 for (int i =0; i<5; ++i) {
  C(i) = complex<double>(i,i);
 }

 for (int i =0; i<2; ++i)
  for (int j=0; j<3; ++j) 
    { 
     A(i,j) = 10*i+ j;
     D (i,j) = A(i,j) / 10.0;
    }

 Af = A;

 cout<<" A= "<<A<<endl;
 cout<<" D= "<<D<<endl;
 cout<<" C= "<<C<<endl;
 cout<<" Arange(0,1),range(1,3)  = "<< A(range(),range(1,3))<<endl;

 H5::H5File file( "ess.h5", H5F_ACC_TRUNC );
 h5_write(file,"A",A);
 h5_write(file,"Af",Af);
 h5_write(file,"C",C);
 h5_write(file,"D",D);

 file.createGroup("G");
 h5_write(file,"G/A",A);

 H5::Group G = file.openGroup("G");
 h5_write(G, "A2",A);

 h5_read (file, "A",B);   cout<< "B = "<< B<<endl;
 h5_read (file, "Af",Bf); cout<< "Bf = "<< Bf<<endl;
 h5_read (file, "D",D2);  cout<< "D = "<< D2<<endl;
 h5_read (file, "C",C2);  cout<< "C = "<< C2<<endl;

 //array<long,1> E; h5_read (file, "A",E);   cout<< "E = "<< E<<endl;

 } 
 catch( const char * err) { cout<<err<<endl;}

 return 0;
}


