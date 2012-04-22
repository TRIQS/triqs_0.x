// PUT AT THE TOP !!!!!
#define BOOST_RESULT_OF_USE_DECLTYPE
#include <boost/utility/result_of.hpp>

#include <triqs/utility/proto/algebra.hpp>
#include <boost/type_traits/is_complex.hpp> 
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/arrays/expressions/matrix_algebra.hpp>
//#include <triqs/arrays/expressions/array_algebra.hpp>
#include <vector>
#include <iostream>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>

namespace proto = boost::proto; namespace mpl = boost::mpl; namespace tup = triqs::utility::proto; 
namespace bf = boost::fusion;
using triqs::arrays::matrix;
using triqs::arrays::range;
using triqs::arrays::array;

#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl;

struct my_matrix_valued_function {
 triqs::arrays::array<double, 3> data;
 my_matrix_valued_function (double x): data(2,2,5) { data()=0; for (size_t u=0; u< 5; ++u) {data(0,0,u) = x*u;data(1,1,u) = -x*u;}}
 my_matrix_valued_function(my_matrix_valued_function const &x): data(x.data) { std::cerr  << "COPY my_matrix_valued_function"<<std::endl ; }
 typedef triqs::arrays::matrix_view<double>  result_type;
 triqs::arrays::matrix_view<double> operator()(size_t i) const { return data(range(),range(),i);}
};

// a trait to identity this type 
template <typename T> struct is_a_m_f                            : mpl::false_{};
template <>           struct is_a_m_f<my_matrix_valued_function> : mpl::true_ {};

// a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
template <typename T> struct is_scalar_or_element : mpl::or_< triqs::arrays::expressions::matrix_algebra::IsMatrix<T>, triqs::utility::proto::is_in_ZRC<T> > {};

struct LeafGrammar   : proto::and_< proto::terminal<proto::_>, proto::if_<is_a_m_f<proto::_value>()> > {}; 
struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<is_scalar_or_element<proto::_value>()> > {}; 

struct myeval {
 BOOST_PROTO_CALLABLE();

 template<typename Sig> struct result;
 template<typename This, typename F, typename AL> struct result<This(F,AL)>
  : boost::result_of<typename boost::remove_reference<F>::type 
  (typename bf::result_of::at< typename boost::remove_reference<AL>::type , mpl::int_<0> >::type )> {};

 template<typename F, typename AL>
  typename boost::result_of<F(typename bf::result_of::at<AL const, mpl::int_<0> >::type)>::type 
  operator ()(F const &f, AL const & al) const { return f(bf::at<mpl::int_<0> >(al)); }
};

struct Grammar : 
 proto::or_<
 proto::when< ScalarGrammar,                       proto::_value>
 ,proto::when< LeafGrammar,                        myeval(proto::_value, proto::_state) >
 ,proto::when< proto::binary_expr <proto::_,Grammar,Grammar>,      proto::_default<Grammar> > 
/* ,proto::when< proto::plus <Grammar,Grammar>,      proto::_default<Grammar> > 
 ,proto::when< proto::minus <Grammar,Grammar>,     proto::_default<Grammar> >
 ,proto::when< proto::multiplies<Grammar,Grammar>, proto::_default<Grammar> >
 ,proto::when< proto::divides<Grammar,Grammar>,    proto::_default<Grammar> >
*/
  ,proto::when< proto::negate<Grammar >,            proto::_default<Grammar> >
 > {};

struct TL : 
 proto::or_<
 proto::when< ScalarGrammar,                       proto::_value>
 ,proto::when< LeafGrammar,                        myeval(proto::_value, proto::_state) >
 ,proto::when< proto::binary_expr <proto::_,TL,TL>,      proto::_default<Grammar> > 
 ,proto::when< proto::unary_expr<proto::_,TL >,               proto::_default<Grammar> >
 > {};




template <typename Expr> struct The_Expr;  // the expression

typedef tup::domain<Grammar,The_Expr,true>  expr_domain; // the domain 

template<typename Expr> struct The_Expr : boost::proto::extends<Expr, The_Expr<Expr>, expr_domain> { // impl the expression
 The_Expr( Expr const & expr = Expr() ) : boost::proto::extends<Expr, The_Expr<Expr>, expr_domain> ( expr ) {}
 template<typename T> typename boost::result_of<Grammar(Expr,bf::vector<T>) >::type 
  operator() (T const & x) const { return Grammar()(*this, bf::make_vector(x)); }
 friend std::ostream &operator <<(std::ostream &sout, The_Expr<Expr> const &expr) { return boost::proto::eval(expr, triqs::utility::proto::algebra::print_ctx (sout)); }
};

BOOST_PROTO_DEFINE_OPERATORS(is_a_m_f,expr_domain);

template<typename T> void print(T const & x) { std::cout << "any T"<<std::endl ;}
template<typename A> void print(The_Expr<A> const & x) { std::cout << "an expression "<< x << std::endl ;}

int main() { 

 my_matrix_valued_function f1(1),f2(10);
 triqs::arrays::matrix<double> A(2,2); A()=0; A(0,0) = 2; A(1,1) = 20;

 TEST( (f1 - f2));

 TEST(f1(1));

 TEST( (2.0* f1 ) (0));
 TEST( (2* f1 + f2) (1));

 TEST( triqs::arrays::matrix<double> ( ( 2*f1  +f2 ) (1)) );
 TEST( triqs::arrays::eval ( ( 2*f1  +f2 ) (1)) );
 TEST( triqs::arrays::eval ( ( 2*f1  - f2 ) (1)) );
 TEST( triqs::arrays::eval ( ( 2*f1  + 8*f2 ) (1)) );
 TEST( triqs::arrays::eval ( ( 2*f1  - f2  - f2) (1)) );
 TEST( triqs::arrays::eval ( ( 2*f1  + 8*f2 + f2 ) (1)) );
 TEST( triqs::arrays::eval ( ( 2*f1  - f2  - 3*f2) (1)) );
 TEST( triqs::arrays::eval ( ( 2*f1  + 8*f2 + f2/2 ) (1)) );
 TEST( triqs::arrays::eval ( ( A*f1  +f2 ) (1)) );

 print ( f1);
 print ( f1 - f2);

};


