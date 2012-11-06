
#include "./common.hpp"

int main() { 

 { 
  BOOST_AUTO (r,  x_ >> 2 );
  //std::cout<<"expr = "<< r << triqs::utility::typeid_name(r)<< std::endl;
 
  BOOST_AUTO (r3,  r(3) );
  std::cout<<"expr = "<< r3 << std::endl;
  std::cout<<"-------------"<<std::endl;
  
  }
{ 
  BOOST_AUTO( expr , 2*x_  + 1);
  BOOST_AUTO (r,  make_function(expr, x_) );
  BOOST_AUTO (r4,  x_ >> expr );
  std::cerr<<"expr = "<< r4 << std::endl;
  //std::cout<<"expr = "<< r << triqs::utility::typeid_name(r)<< std::endl;
 
  BOOST_AUTO (r3,  r(3) );
  std::cout<<"expr = "<< r3 << std::endl;
  std::cout<<"-------------"<<std::endl;
  
  }
 { 

  TEST( (y_ >>  x_  + 2*y_) );
  TEST( x_ >> (y_ >>  x_  + 2*y_ ));
  //TEST(( y_ >> x_ >>  x_  + 2*y_ ));

  BOOST_AUTO( r,y_ >>  x_  + 2*y_ );
   std::cout<<" r >> = "<< r << std::endl;
  //std::cout << triqs::utility::typeid_name(r2)<< std::endl;
  // std::cout<<" r >> = "<< r << triqs::utility::typeid_name(r)<< std::endl;
 
  BOOST_AUTO( r2, x_ >> r );
  std::cout<<" r2 >> = "<< r2 << std::endl;
  //std::cout << triqs::utility::typeid_name(r2)<< std::endl;
  //std::cout<<"-------------"<<std::endl;
  
  //BOOST_AUTO( r3, x_ >> y_ >>  x_  + 2*y_  );
 } 
 { 
  BOOST_AUTO( expr , x_  + 2*y_);

  std::cout<<"expr = "<< expr<< std::endl;
  std::cout << " ph list of expr "<< placeholder_list_as_string(expr)<<std::endl;

  BOOST_AUTO (r,  make_function(expr, x_) );
  
  //std::cout<<"expr = "<< r << triqs::utility::typeid_name(r)<< std::endl;
 
  std::cout << " ph list of r "<< placeholder_list_as_string(r)<<std::endl;

  BOOST_AUTO (r2,  tql::eval( r , y_ = 2) );
  //std::cout<<"expr = "<< r2 << triqs::utility::typeid_name(r2)<< std::endl;
  std::cout << " ph list of r2 "<< placeholder_list_as_string(r2)<<std::endl;
  
  BOOST_AUTO (r3,  r2(3) );
  std::cout<<"expr = "<< r3 << std::endl;
  
  TEST( tql::eval ( make_function(x_ + 2*y_ , x_) , y_ = 2) (3));

  TEST( make_function(x_ + 2 , x_)  (3));

  std::cout<<"-------------"<<std::endl;
 }

{ 
  BOOST_AUTO( expr , x_  + 2*y_ + z_);

  std::cout<<"expr = "<< expr<< std::endl;
  std::cout << " ph list of expr "<< placeholder_list_as_string(expr)<<std::endl;

  BOOST_AUTO (r,  make_function(expr, x_) );
  std::cout << " ph list of r "<< placeholder_list_as_string(r)<<std::endl;

  BOOST_AUTO (r2,  tql::eval(r , y_ = 2) );
  std::cout << " ph list of r2 "<< placeholder_list_as_string(r2)<<std::endl;
  
  std::cout<<"-------------"<<std::endl;
 }



}

