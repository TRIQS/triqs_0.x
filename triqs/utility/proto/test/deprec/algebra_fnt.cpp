#include <triqs/utility/proto/algebra.hpp>
#include <triqs/utility/proto/tools.hpp>
#include <boost/type_traits/is_complex.hpp> 
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/arrays/expressions/matrix_algebra.hpp>
//#include <triqs/arrays/expressions/array_algebra.hpp>
#include <vector>
#include <iostream>

namespace tqa=triqs::arrays; namespace mpl = boost::mpl; namespace tup = triqs::utility::proto; 
using tqa::matrix;
using tqa::range;
using tqa::array;

#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl;

struct my_matrix_valued_function {
 tqa::array<double, 3> data;
 my_matrix_valued_function (double x): data(2,2,5) { data()=0; for (size_t u=0; u< 5; ++u) {data(0,0,u) = x*u;data(1,1,u) = -x*u;}}
 my_matrix_valued_function(my_matrix_valued_function const &x): data(x.data) { std::cerr  << "COPY my_matrix_valued_function"<<std::endl ; }
 tqa::matrix_view<double> operator()(size_t i) const { return data(range(),range(),i);}
};

// a trait to identity this type 
template <typename T> struct is_a_m_f                            : mpl::false_{};
template <>           struct is_a_m_f<my_matrix_valued_function> : mpl::true_ {};

// a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
template <typename T> struct is_scalar_or_element : mpl::or_< tqa::is_matrix_expr<T>, triqs::utility::proto::is_in_ZRC<T> > {};

// This macro declares the algebra of algebra-valued functions...
//TRIQS_PROTO_DEFINE_ALGEBRA_VALUED_FNT_ALG (is_a_m_f, is_scalar_or_element);
// second solution is better since I can reuse T for pattern matching...
typedef tup::domain_and_expression_generator< tup::algebra::grammar_generator1< tup::algebra::algebra_function_desc,is_a_m_f, is_scalar_or_element >::type  > _D;
BOOST_PROTO_DEFINE_OPERATORS(is_a_m_f, _D::expr_domain);

template<typename T> void print(T const & x) { std::cout << "any T"<<std::endl ;}
template<typename A> void print(_D::The_Expr<A> const & x) { std::cout << "an expression "<< x << std::endl ;}

int main() { 

 my_matrix_valued_function f1(1),f2(10);
 tqa::matrix<double> A(2,2); A()=0; A(0,0) = 2; A(1,1) = 20;

 TEST( (f1 - f2));

 TEST(f1(1));

 TEST( (2.0* f1 ) (0));
 TEST( (2* f1 + f2) (1));

 TEST( tqa::make_matrix ( ( 2*f1  +f2 ) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  +f2 ) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  - f2 ) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  + 8*f2 ) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  - f2  - f2) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  + 8*f2 + f2 ) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  - f2  - 3*f2) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  + 8*f2 + f2/2 ) (1)) );
 TEST( tqa::make_matrix ( ( A*f1  +f2 ) (1)) );
 print ( f1);
 print ( f1 - f2);
};


