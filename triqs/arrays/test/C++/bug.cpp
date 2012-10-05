#include "./python_stuff.hpp"
#include "./src/array.hpp"
#include "./src/matrix.hpp"
#include "./src/linalg/matmul.hpp"
#include <iostream>
#include <triqs/arrays/proto/array_algebra.hpp>
#include <triqs/arrays/proto/matrix_algebra.hpp>

using std::cout; using std::endl;
using namespace triqs::arrays;
using namespace triqs::arrays;

int main() { 
  typedef Option::C storage_order;
  // typedef Option::Fortran storage_order;
 
 array<std::complex<double>,3,storage_order> TMPALL (2,2,5); TMPALL()=-1;
 matrix_view<std::complex<double>,storage_order > TMP ( TMPALL (range(), range(), 2));
 matrix<std::complex<double>,storage_order > M1(2,2); M1()=0; M1(0,0) = 2; M1(1,1) = 3.2;
 TMP() =0; 
 auto R = M1*M1; // +1.0; // * M1; 
 std::cout << R << std::endl ;
 std::cout << R[mini_vector<size_t,2>(0,0)] << R[mini_vector<size_t,2>(1,1)] << std::endl ;
 TMP() =  R; // M1*( M1 ); 
 std::cerr  << "DEBUG TAIL INVERIN TMP " <<TMP<< std::endl; 
 std::cerr  << "DEBUG TAIL INVERIN TMP " <<TMPALL<< std::endl; 

matrix<std::complex<double>,storage_order > TMP2(2,2); 
  TMP2() = R;
 std::cerr  << "DEBUG TAIL INVERIN TM2P " <<TMP2<< std::endl; 

}
