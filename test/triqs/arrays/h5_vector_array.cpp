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
#include "./common.hpp"
#include <iostream>
#include "./src/asserts.hpp"
#include <triqs/h5.hpp>

using std::cout; using std::endl;
namespace tqa = triqs::arrays;
namespace h5 = triqs::h5;
using tqa::range;
using namespace tqa;

int main(int argc, char **argv) {

#ifndef TRIQS_WORKAROUND_INTEL_COMPILER_BUGS
 try { 

  std::vector <double> v {1.1,2.2,3.3,4.5};
  std::vector <std::complex<double>> vc {1.1,2.2,3.3,4.5};

  std::vector <double> v2;
  std::vector <std::complex<double>> vc2;

  {
   H5::H5File file1( "test_std_vector.h5", H5F_ACC_TRUNC );
   h5::group top(file1);
   h5_write(top,"vdouble",v);
   h5_write(top,"vcomplex",vc);
  }

  H5::H5File file2( "test_std_vector.h5", H5F_ACC_RDONLY );
  h5::group top2(file2);

  h5_read(top2,"vdouble",v2);
  h5_read(top2,"vcomplex",vc2);

  for (int i = 0; i <v.size(); ++i) assert_close(v[i],v2[i]); 
  for (int i = 0; i <vc.size(); ++i) assert_close(vc[i],vc2[i]); 
 }
 catch(std::exception const &e) { std::cerr << e.what() << std::endl ;}
#endif

 return 0;
}

