#include "../core.hpp"
#include <iostream> 
#include "triqs/utility/typeid_name.hpp"
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>

#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl;

struct F1 : boost::noncopyable { 
 double v; 
 F1(double v_): v(v_){} 
 double operator() (double x) const { return 10*x;}
 TRIQS_CLEF_ADD_LAZY_CALL_REF(1,F1);

 TRIQS_CLEF_HAS_AUTO_ASSIGN(); 
 template<typename Fnt> friend void triqs_nvl_auto_assign (F1 & x, Fnt f) { x.v++; std::cerr<< " called triqs_nvl_auto_assign "<< f(100)<<std::endl;}
 friend std::ostream & triqs_nvl_formal_print(std::ostream & out, F1 const & x) { return out<<"F1";}

};

struct F1b {
 double v; boost::shared_ptr<double> v_ptr;

 F1b(double v_): v_ptr (new double(v_)) {v=v_;} // on purpose, I avoid adding a default constructor !  
 F1b (F1b const & x)  { v_ptr  = x.v_ptr; v = x.v; std::cerr<< " F1 COPY "<<std::endl;}

 double operator() (double x) const { return 10*x;}
 TRIQS_CLEF_ADD_LAZY_CALL_WITH_COPY(1,F1b);

 TRIQS_CLEF_HAS_AUTO_ASSIGN(); 
 template<typename Fnt> friend void triqs_nvl_auto_assign (F1b & x, Fnt f) { x.v++;(*(x.v_ptr))++; std::cerr<< " called triqs_nvl_auto_assign "<< f(100)<<std::endl;}
 friend std::ostream & triqs_nvl_formal_print(std::ostream & out, F1b const & x) { return out<<"F1b";}

};

struct F2 {

 double v; 
 F2() {v=0;}

 double operator()(double x, double y) const { return 10*x + y;}
 TRIQS_CLEF_ADD_LAZY_CALL_REF(2,F2);

 TRIQS_CLEF_HAS_AUTO_ASSIGN(); 
 template<typename Fnt> friend void triqs_nvl_auto_assign (F2 & x, Fnt f) { x.v++; std::cerr<< " called set_from_function "<< f(10,20)<<std::endl;}
 friend std::ostream & triqs_nvl_formal_print(std::ostream & out, F2 const & x) { return out<<"F2";}
};


using namespace triqs::clef;

triqs::clef::placeholder <1> x_; 
triqs::clef::placeholder <2> y_; 
triqs::clef::placeholder <3> z_; 
namespace tql= triqs::clef;

