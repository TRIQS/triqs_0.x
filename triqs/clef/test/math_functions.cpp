#include "./common.hpp"
#include <boost/function.hpp>
#include "triqs/clef/math.hpp"


double foo(double x) { return x/2;}
int foo(int x) { return x/2;}

double bar(double x, double y) { return x+y;}

namespace triqs { namespace clef { 

 using ::foo; 
 using ::bar; 

 template<typename T> 
 typename boost::disable_if<triqs::clef::is_lazy<T>,T>::type 
 inc (T const & x) { return x+1;}
// using ::inc; 
 // moving the declaration up and using using:: does not work on gcc4.4
 // however not using using:: does not work on icc 11 !

 TRIQS_CLEF_MAKE_FNT_LAZY (2, bar) ;
 TRIQS_CLEF_MAKE_FNT_LAZY (1, inc) ;
 TRIQS_CLEF_MAKE_FNT_LAZY (1, foo) ;

 TRIQS_CLEF_MAKE_FNT_LAZY_NOWRAP (1,cos1, boost::function<double (double) >, boost::function<double (double) > ( static_cast<double(*)(double)>(std::cos)) );
 TRIQS_CLEF_MAKE_FNT_LAZY_NOWRAP (1,cos2, double (*)(double) , static_cast<double(*)(double)>(std::cos)) ;
 TRIQS_CLEF_MAKE_FNT_LAZY_NOWRAP (1,foo1, double (*)(double) , static_cast<double(*)(double)>(::foo)) ;
 TRIQS_CLEF_MAKE_FNT_LAZY_NOWRAP (2,bar1, double (*)(double,double) , static_cast<double(*)(double,double)>(::bar)) ;

}}

namespace tql= triqs::clef;
int main() { 
 using std::cout; using  std::endl;

 TEST  ( tql::eval ( cos(x_) ,x_=2) );
 TEST  ( tql::eval ( cos(2*x_+1) ,x_=2) );
 TEST  ( tql::eval ( cos1(2*x_+1) ,x_=2) );
 TEST  ( tql::eval ( cos2(2*x_+1) ,x_=2) );

 TEST  ( tql::eval ( abs(2*x_-1) ,x_=2) );
 TEST  ( tql::eval ( abs(2*x_-1) ,x_=-2) );
 TEST  ( tql::eval ( floor(2*x_-1) ,x_=2.3) );
 TEST  ( tql::eval ( pow(2*x_+1,2) ,x_=2.0) );

 TEST  ( tql::eval ( foo1(2*x_+1) , x_ = 2) ); 
 TEST  ( tql::eval ( foo(2*x_+1) , x_ = 2) ); 
 TEST  ( tql::eval ( foo(2*x_+1) , x_ = 2.0) ); 
 TEST  ( tql::eval ( bar1(2*x_+1,x_ -1) , x_ = 2) ); 
 TEST  ( tql::eval ( bar(2*x_+1,x_ -1) , x_ = 2) ); 
 TEST  ( tql::eval ( inc(2*x_+1) , x_ = 2) ); 

}

