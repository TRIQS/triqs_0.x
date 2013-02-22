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
 P["s"] = std::string("-14.3");
 P["sc"] = "-14.3";

 long j = P["a"];
 double x = P["d"];
 double y = P["a"];
 double z = P["s"];
 double zc = P["sc"];

 std::cout  << j << std::endl ;
 std::cout  << x << std::endl;
 std::cout  << y << std::endl ; 
 std::cout  << z << std::endl ; 
 std::cout  << zc << std::endl ; 
 std::cout  << P["a"] << std::endl ; 

 // testing that copy is a copy 
 parameters P2 = P; 
 P2["a"] = 12.3;

 // Put P2 in P ... 
 P["P2"] = P2;

 std::cout  << P << std::endl ;
 {
  H5::H5File file( "ess.h5", H5F_ACC_TRUNC );
  h5_write( file, "Parameters", P);
 }

 {
  H5::H5File file( "ess2.h5", H5F_ACC_TRUNC );
  h5_write( file, "Parameters", P2);
 }


 std::string s = triqs::serialize(P);
 //std::cout  << " serialization "<< s << std::endl; 

 parameters P3 = triqs::deserialize<parameters>(s);
 {
  H5::H5File file( "ess3.h5", H5F_ACC_TRUNC );
  h5_write( file, "Parameters", P3);
 }

 parameters P4;
 std::cout << "P4 before : "<< P4<< std::endl ;
 {
  H5::H5File file( "ess2.h5", H5F_ACC_RDONLY );
  h5_read( file, "Parameters", P4);
 }
 {
  H5::H5File file( "ess2_relo.h5", H5F_ACC_TRUNC );
  h5_write( file, "Parameters", P4);
 }
std::cout << "P4 after : "<< P4<< std::endl ;

 return 0;
}
