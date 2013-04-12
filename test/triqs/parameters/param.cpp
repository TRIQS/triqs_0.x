#include <triqs/parameters.hpp>
#include <iostream>
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

using namespace triqs;
using namespace triqs::utility;

// what is the concept of things that can be put in the dict ?
struct my_obj { 
 int i;
 my_obj(int i_) { i=i_;}

 my_obj() { i=0;}
 my_obj(my_obj const & x) { i = x.i;}
 my_obj(my_obj && x) { std::swap(i,x.i);}

 // h5 operation 
 friend std::string get_triqs_hdf5_data_scheme( my_obj const&) { return "myobject_is_nice";} 

 friend void h5_write ( h5::group F, std::string const & subgroup_name, my_obj const & p){
  h5::group gr =  F.create_group(subgroup_name);
  gr.write_triqs_hdf5_data_scheme(p);
  h5_write(gr,"data",p.i);
 }

 friend void h5_read ( h5::group F, std::string const & subgroup_name, my_obj & p){
  auto gr = F.open_group(subgroup_name);
  h5_read(gr,"data",p.i);
 }

 // print
 friend std::ostream & operator << (std::ostream & out, my_obj const & p) { return out<< "myobject("<<p.i<<")";}

 // boost serialization 
 friend class boost::serialization::access;
 template<class Archive>
  void serialize(Archive & ar, const unsigned int version) { ar & boost::serialization::make_nvp("i",i); }
};


int main() {

 parameters P;

 P["myobject1"] = my_obj(18);

 P["a"] = long(1);
 P["d"] = 2.7;
 P["s"] = std::string("-14.3");
 P["sc"] = "-14.3";

 triqs::arrays::array<double,2> A(2,2); A()=0;A(0,0) = 1.3; A(1,1) = -8.2;
 triqs::arrays::array<long,1> B(3); B()=0;B(0) = 3; B(1) = -8;
 P["A"] = std::move(A);
 P["B"] = B;
 std::cout  << "A"<< P["A"] << std::endl;
 std::cout  << "B"<< P["B"] << std::endl;

 triqs::arrays::array<long,1> C;
 C = extract<decltype(C)>(P["B"]);
 std::cout  << "C" << C << std::endl;

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

 std::cout  << P << std::endl;

 {
  H5::H5File file( "ess.h5", H5F_ACC_TRUNC );
  h5_write( file, "Parameters", P);
 }

 {
  H5::H5File file( "ess2.h5", H5F_ACC_TRUNC );
  h5_write( file, "Parameters", P2);
 }

 std::string s = triqs::serialize(P);
 parameters P3 = triqs::deserialize<parameters>(s);
 {
  H5::H5File file( "ess3.h5", H5F_ACC_TRUNC );
  h5_write( file, "Parameters", P3);
 }

 parameters P4;
 std::cout << "P4 before : "<< P4<< std::endl ;
 {
  H5::H5File file( "ess.h5", H5F_ACC_RDONLY );
  h5_read( file.openGroup("/"), "Parameters", P4);
 }
 {
  H5::H5File file( "ess_relo.h5", H5F_ACC_TRUNC );
  h5_write( file.openGroup("/"), "Parameters", P4);
 }
 std::cout << "P4 after : "<< P4<< std::endl ;

 return 0;
}
