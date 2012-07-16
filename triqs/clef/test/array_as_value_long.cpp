#include "./common.hpp"
#include "triqs/clef/adapters/array.hpp"

triqs::clef::placeholder <4> i_; 
triqs::clef::placeholder <5> j_; 
namespace tqa= triqs::arrays;

int main() { 

 //double x=1,y=2;

 triqs::arrays::array<double,1> A(3); A(0)=1; A(1)=2; A(2)=3;
 triqs::arrays::array_view<double,1> V(A); 
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
 B(j_,i_) = 10*i_ + j_ +1; 
 B(j_,i_) = 10*i_ + 2*j_; 
 B(j_,i_) = 10*i_ - j_; 
 B(j_,i_) = i_ + j_; 
 B(j_,i_) = 1+2+ 10*i_ + j_; 
 B(j_,i_) = 10*i_ + j_;
 B(j_,i_) = 10*i_ + j_ +1;
 B(j_,i_) = 10*i_ + 2*j_;
 B(j_,i_) = 10*i_ - j_;
 B(j_,i_) = i_ + j_;   
 B(j_,i_) = 1+2+ 10*i_ + j_;
 B(j_,i_) = 10*i_ + j_/1.0;
  B(j_,i_) = 10*i_ + j_/3.8 +1;
   B(j_,i_) = 10*i_ + 2*j_/8.1;
    B(j_,i_) = 10*i_ - j_/9.0;
     B(j_,i_) = i_ + j_/8.1;
      B(j_,i_) = 1+2+ 10*i_ + j_/7.2;

 B(j_,i_) = 1+2+ 10*i_ + j_/2.0;
 B(j_,i_) = 11*i_ + j_; 
 B(j_,i_) = 12*i_ + j_; 
 B(j_,i_) = 13*i_ + j_; 
 B(j_,i_) = 14*i_ + j_; 
 B(j_,i_) = 15*i_ + j_; 
 B(j_,i_) = 20*i_ + j_; 
 B(j_,i_) = 30*i_ + j_; 
 B(j_,i_) = 40*i_ + j_; 
 B(j_,i_) = 50*i_ + j_; 
 std::cout <<B<<std::endl;


}

