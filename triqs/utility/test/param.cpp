#include <triqs/utility/parameters.hpp>
#include <iostream>

using namespace triqs::utility;

int main() {

 parameters P;

/* P 
  ( "A", long(0), " comment ... ")
  ( "A", long(0), " comment ... ")
  ( "A", long(0), " comment ... ")
  ( "A", long(0), " comment ... ")
  ( "A", long(0), " comment ... ")
  ( "A", long(0), " comment ... ")
  ;
*/

 P["a"] = long(1);
 P["d"] = 2.7;
 
 long j = P["a"];
 double x = P["d"];
 double y = P["a"];
 	
 std::cout  << x << std::endl;
 std::cout  << y << std::endl ; 
 std::cout  << j << std::endl ;

 
 H5::H5File file( "ess.h5", H5F_ACC_TRUNC );
 h5_write( file, "a", P.object_map["a"]);
 h5_write( file, "d", P.object_map["d"]);

 return 0;
}
