#include <triqs/utility/proto/algebra.hpp>
#include <boost/type_traits/is_complex.hpp> 
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/arrays/expressions/matrix_algebra.hpp>
//#include <triqs/arrays/expressions/array_algebra.hpp>
#include <vector>
#include <iostream>

namespace mpl = boost::mpl;
using triqs::arrays::matrix;
using triqs::arrays::range;
using triqs::arrays::array;

#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl;

struct my_matrix_valued_function {
 triqs::arrays::array<double, 3> data;
 my_matrix_valued_function (double x): data(2,2,5) { data()=0; for (size_t u=0; u< 5; ++u) {data(0,0,u) = x*u;data(1,1,u) = -x*u;}}
 my_matrix_valued_function(my_matrix_valued_function const &x): data(x.data) { std::cerr  << "COPY my_matrix_valued_function"<<std::endl ; }
 triqs::arrays::matrix_view<double> operator()(size_t i) const { return data(range(),range(),i);}
};

// a trait to identity this type 
template <typename T> struct is_a_m_f                            : mpl::false_{};
template <>           struct is_a_m_f<my_matrix_valued_function> : mpl::true_ {};

// a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
template <typename T> struct is_scalar_or_element   : mpl::or_< triqs::arrays::expressions::matrix_algebra::IsMatrix<T>, triqs::utility::proto::is_in_ZRC<T> > {};

// This macro declare the algebra of algebra-valued functions...
TRIQS_PROTO_DEFINE_ALGEBRA_VALUED_FNT_ALG (is_a_m_f, is_scalar_or_element);

int main() { 

 my_matrix_valued_function f1(1),f2(10);
 triqs::arrays::matrix<double> A(2,2); A()=0; A(0,0) = 2; A(1,1) = 20;

 TEST( (f1 - f2));

 TEST(f1(1));

 TEST( (2.0* f1 ) (0));
 TEST( (2* f1 + f2) (1));

 TEST( triqs::arrays::matrix<double> ( ( 2*f1  +f2 ) (1)) );
 TEST( triqs::arrays::eval ( ( 2*f1  +f2 ) (1)) );

 TEST( triqs::arrays::eval ( ( A*f1  +f2 ) (1)) );

};


