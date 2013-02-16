#include "./common.hpp"
#include "triqs/clef/adapters/array.hpp"

triqs::clef::placeholder <4> i_; 
triqs::clef::placeholder <5> j_; 
namespace tqa= triqs::arrays;

int main() { 

 //double x=1,y=2;

 tqa::array<double,1> A(3); A(0)=1; A(1)=2; A(2)=3;
 tqa::array_view<double,1> V(A); 
 //double a = 2;

 BOOST_AUTO( exp , x_ + A);

 std::cerr<< " -------------"<<std::endl;
 TEST( tqa::make_array( tql::eval( 2.0*x_ , x_ =A) ));
 TEST( tqa::make_array( tql::eval(2.0*x_ , x_ =3*A) ));


 TEST( tqa::make_array( tql::eval (A + 2*x_ , x_ =A) )); 
 TEST( ( tql::eval (A + 2*x_ , x_ =A) )); 

 TEST( A(i_));
 TEST( tql::eval( A(i_), i_=1));

 triqs::arrays::array<double,2> B(2,2);
 B(1,1) = 10;
 TEST( tql::eval( B(i_,j_) , i_ = 1, j_ = 1)); 

 B(i_,j_) = 10*i_ + j_; 
 std::cout <<B<<std::endl;

 B(j_,i_) = 10*i_ + j_; 
 std::cout <<B<<std::endl;

 B(i_,j_) = i_ - 2* j_; 
 std::cout <<B<< std::endl;

 //TEST( A + 2*A);

}

