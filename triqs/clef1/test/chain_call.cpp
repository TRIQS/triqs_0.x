#include "./common.hpp"

#include <vector>

struct F2_vec { 
 std::vector<F2> vec; 
 F2_vec(size_t siz=0): vec (siz) {}
 F2_vec(F2_vec const & x): vec(x.vec) {}
 
 F2 const & operator()(size_t i) const { return vec[i];}
 TRIQS_CLEF_ADD_LAZY_CALL_WITH_COPY(1,F2_vec);

 F2 const & operator[](size_t i) const { return vec[i];}
 TRIQS_CLEF_ADD_LAZY_SUBSCRIPT_WITH_COPY(F2_vec);

 TRIQS_CLEF_HAS_AUTO_ASSIGN_SUBSCRIPT(); 
 template<typename Fnt>
  friend void triqs_nvl_auto_assign_subscript (F2_vec & x, Fnt f) { 
   for (size_t i=0; i< x.vec.size(); ++i) triqs_nvl_auto_assign(x.vec[i],  f(i)); } 

 friend std::ostream & triqs_nvl_formal_print(std::ostream & out, F2_vec const & x) { return out<<"F2_vec";}
};

struct F2_vecB { 
 std::vector<F2> vec; 
 F2_vecB(size_t siz=0): vec (siz) {}
 F2_vecB(F2_vecB const & x): vec(x.vec) {}
 
 F2 const & operator()(size_t i) const { return vec[i];}
 TRIQS_CLEF_ADD_LAZY_CALL_WITH_COPY(1,F2_vecB);

 F2 const & operator[](size_t i) const { return vec[i];}
 TRIQS_CLEF_ADD_LAZY_SUBSCRIPT_WITH_COPY(F2_vecB);

 TRIQS_CLEF_HAS_AUTO_ASSIGN(); 
 template<typename Fnt>
  friend void triqs_nvl_auto_assign (F2_vecB & x, Fnt f) { 
   for (size_t i=0; i< x.vec.size(); ++i) triqs_nvl_auto_assign(x.vec[i],  f(i)); } 

 friend std::ostream & triqs_nvl_formal_print(std::ostream & out, F2_vecB const & x) { return out<<"F2_vec";}
};

triqs::clef::placeholder <4> i_; 

using  namespace triqs::clef;
int main() { 

 F2_vec V(3);
 F2_vecB Vb(3);

 //double x=1,y=2;

 { 
  BOOST_AUTO( expr , Vb(i_)(x_,y_) );
  std::cout<<"expr = "<< expr<< std::endl;

  TEST( tql::eval(expr, x_=2,y_=3, i_=0));
  TEST( tql::eval(expr, x_=2));
  TEST( tql::eval(expr, x_=2,i_=1));
  TEST( tql::eval(expr, x_=2,y_=3));
  TEST( tql::eval(  tql::eval(expr, x_=2,y_=3) , i_=0) );
  std::cout<<"-------------"<<std::endl;
 }

 {
  BOOST_AUTO( expr , V[i_](x_,y_) );
  std::cout<<"expr = "<< expr<< std::endl;

  TEST( tql::eval(expr, x_=2,y_=3, i_=0));
  TEST( tql::eval(expr, x_=2));
  TEST( tql::eval(expr, x_=2,i_=1));
  TEST( tql::eval(expr, x_=2,y_=3));
  TEST( tql::eval(  tql::eval(expr, x_=2,y_=3) , i_=0) );
  std::cout<<"-------------"<<std::endl;

 }
 {
  // test assign
  Vb(i_)(x_,y_)  = x_ + 10*y_ + 100*i_;
 }

 {
  // test assign
  V[i_](x_,y_)  = x_ + 10*y_ + 100*i_;
 }

}

