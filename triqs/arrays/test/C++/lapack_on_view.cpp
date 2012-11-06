#include "./python_stuff.hpp"
#include "./src/array.hpp"
#include "./src/matrix.hpp"
#include "./src/linalg/matmul.hpp"
#include <iostream>
#include <triqs/arrays/proto/array_algebra.hpp>
#include <triqs/arrays/proto/matrix_algebra.hpp>
#include <triqs/arrays/asserts.hpp>

using std::cout; using std::endl;
using namespace triqs::arrays;
using namespace triqs::arrays;

// to be extended to more complex case
// calling lapack on view to test cache securities....
//
int main() { 
  typedef Option::C storage_order;
  // typedef Option::Fortran storage_order;
 
 array<std::complex<double>,3,storage_order> TMPALL (2,2,5); TMPALL()=-1;
 matrix_view<std::complex<double>,storage_order > TMP ( TMPALL (range(), range(), 2));
 matrix<std::complex<double>,storage_order > M1(2,2), Res(2,2); 
 M1()=0; M1(0,0) = 2; M1(1,1) = 3.2;
 Res()=0; Res(0,0) = 8; Res(1,1) = 16.64;
 TMP() =0; 
 TMP() =   M1*( M1 + 2.0 ); 
 assert_all_close(TMP(), Res, 1.e-10); 
}
