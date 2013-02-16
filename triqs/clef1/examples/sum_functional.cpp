#include <triqs/clef/core.hpp>
#include <iostream>
namespace tql = triqs::clef;

// a general implementation type for all domain and summation methods
template<typename Domain, typename Options> struct sum_impl;

struct Riemann{}; // one sum method

template<typename Domain> struct sum_impl<Domain, Riemann> {
 Domain d;
 sum_impl(Domain const &d_):d(d_) {}

 // Disabled when F is a lazy_expression
 template <typename F>
  typename boost::disable_if< tql::is_lazy <F>, double >::type 
  operator() (F const & f) const { double s=0; for (int u=0; u<10; ++u) s += f(u/10.0); return s;}

 // Defines operator(). Enabled iff one argument is a clef expression.
 TRIQS_CLEF_ADD_LAZY_CALL_WITH_COPY(1,sum_impl);

 // How to print this object in the printing of clef expression
 friend std::ostream & triqs_nvl_formal_print(std::ostream & out, sum_impl const & x) { return out<<"sum";}
};

template<typename Domain, typename Option> 
sum_impl<Domain,Option> sum_functional (Domain const & d, Option) { return sum_impl<Domain,Option> (d);}

//---------    MAIN   ---------------------------------------

struct DOM{};

int main() { 

 // two placeholders
 tql::placeholder <1> x_; 
 tql::placeholder <2> y_; 

 DOM d; // a domain
 
 // integrate_on_d is the integration (!) over d with Riemann method
 BOOST_AUTO( integrate_on_d , sum_functional( d, Riemann()) ); 
 
 // This is a simple application of the sum to a function 
 std::cout<< integrate_on_d( x_ >> 2*x_ + 1 )  << std::endl;

 // This creates a clef expression of placeholder y_, waiting to make the sum
 // Indeed the argument is a clef expression of y_ returning a function
 // composed by integrate_on_d, it is a clef expression of y_ returning a double
 std::cout<< integrate_on_d( x_ >> 2*x_ + 1 + y_ )  << std::endl;

 // Of course this expression can be mixed with others...
 std::cout<< y_ + 2* integrate_on_d( x_ >> 2*x_ + 1 + y_ )  << std::endl;

 // and it can be evaluated 
 std::cout<< tql::eval (y_ + 2* integrate_on_d( x_ >> 2*x_ + 1 + y_ ) ,y_ =0)  << std::endl;
 
}

