#include <triqs/utility/proto/tools.hpp>
#include <boost/type_traits/is_complex.hpp> 
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/arrays/expressions/matrix_algebra.hpp>
//#include <triqs/arrays/expressions/array_algebra.hpp>
#include <vector>
#include <iostream>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>

namespace tqa=triqs::arrays; 
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
 typedef bool D;
 bool b() const { return true;}
 friend std::ostream & formal_print(std::ostream & out, my_matrix_valued_function const & x) { return out<<"my_matrix_valued_function";}
 
};

// a trait to identity this type 
template <typename T> struct is_a_m_f                            : mpl::false_{};
template <>           struct is_a_m_f<my_matrix_valued_function> : mpl::true_ {};

// a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
template <typename T> struct is_scalar_or_element : mpl::or_< triqs::arrays::ImmutableMatrix<T>, tup::is_in_ZRC<T> > {};

struct ElementGrammar: proto::and_< proto::terminal<proto::_>, proto::if_<is_a_m_f<proto::_value>()> > {}; 
struct ScalarGrammar : proto::and_< proto::terminal<proto::_>, proto::if_<is_scalar_or_element<proto::_value>()> > {}; 

struct algebra_grammar : 
 proto::or_<
 ScalarGrammar , ElementGrammar
 , proto::plus      <algebra_grammar,algebra_grammar>
 , proto::minus     <algebra_grammar,algebra_grammar>
 , proto::multiplies<algebra_grammar,algebra_grammar>
 , proto::divides   <algebra_grammar,algebra_grammar>
 , proto::negate    <algebra_grammar >
 > {};

template<int Arity> struct eval_fnt;
template<> struct eval_fnt<1> {
 BOOST_PROTO_CALLABLE();
 template<typename Sig> struct result;
 template<typename This, typename F, typename AL> struct result<This(F,AL)>
  : boost::result_of<typename boost::remove_reference<F>::type 
  (typename bf::result_of::at< typename boost::remove_reference<AL>::type , mpl::int_<0> >::type )> {};

 template<typename F, typename AL>
  typename boost::result_of<F(typename bf::result_of::at<AL const, mpl::int_<0> >::type)>::type 
  operator ()(F const &f, AL const & al) const { return f(bf::at<mpl::int_<0> >(al)); }
};

template<typename Evaluator> struct eval_transform : proto::or_<
 proto::when<ScalarGrammar, proto::_value>
 ,proto::when<ElementGrammar, Evaluator(proto::_value, proto::_state) >
 ,proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::_default<eval_transform<Evaluator> > >
 //,proto::when< proto::binary_expr <proto::_,eval_transform<Evaluator>,eval_transform<Evaluator> >, proto::_default<eval_transform<Evaluator> > > 
 //,proto::when< proto::negate<eval_transform<Evaluator> >, proto::_default<eval_transform<Evaluator> > >
 > {};

template<typename Accumulator> struct accumulator_transform : proto::or_<
 proto::when< ScalarGrammar, proto::_state>
 ,proto::when< ElementGrammar, Accumulator(proto::_value, proto::_state)  >
 , proto::when<proto::nary_expr<proto::_, proto::vararg<proto::_> >,  proto::fold<proto::_, proto::_state, accumulator_transform<Accumulator> >() >
 // ,proto::when< proto::binary_expr <proto::_,accumulator_transform,accumulator_transform>, accumulator_transform(proto::_left, accumulator_transform( proto::_right)) >
 // ,proto::when< proto::unary_expr<proto::_,accumulator_transform >,     accumulator_transform(proto::_left) >
 > {};

//---------------------------

template <typename Expr> struct expr_a_m_f;  // the expression

typedef tup::domain<algebra_grammar,expr_a_m_f,true>  domain_a_m_f; // the domain 

template<typename Expr> struct expr_a_m_f : boost::proto::extends<Expr, expr_a_m_f<Expr>, domain_a_m_f> { // impl the expression
 expr_a_m_f( Expr const & expr = Expr() ) : boost::proto::extends<Expr, expr_a_m_f<Expr>, domain_a_m_f> ( expr ) {}
 template<typename T> typename boost::result_of<eval_transform<eval_fnt<1> >(Expr,bf::vector<T>) >::type 
  operator() (T const & x) const { return eval_transform<eval_fnt<1> >()(*this, bf::make_vector(x)); }
 friend std::ostream &operator <<(std::ostream &sout, expr_a_m_f<Expr> const &expr) { return boost::proto::eval(expr, tup::algebra_print_ctx (sout)); }
};

BOOST_PROTO_DEFINE_OPERATORS(is_a_m_f,domain_a_m_f);

template<typename T> void print(T const & x) { std::cout << "any T"<<std::endl ;}
template<typename A> void print(expr_a_m_f<A> const & x) { std::cout << "an expression "<< x << std::endl ;}

struct myaccu1{
 BOOST_PROTO_CALLABLE();
 template<typename Sig> struct result;
 template<typename This, typename A, typename BFL> struct result<This(A,BFL)> { 
  typedef bf::cons<typename boost::remove_reference<A>,typename boost::remove_reference<BFL>::type> type;
 };
 template<typename A, typename BFL> bf::cons<A,BFL> operator ()(A &a, BFL const & bfl) const { return bf::cons<A,BFL>(a,bfl); }
 //template<typename A, typename BFL> bf::cons<A,BFL> operator ()(A &a, BFL const & bfl BOOST_PROTO_DISABLE_IF_IS_CONST(Seq)) const { return bf::cons<A,BFL>(a,bfl); }
};

struct myaccu{
 BOOST_PROTO_CALLABLE();	
 template<typename Sig> struct result;
 template<typename This, typename A, typename BFL> 
  struct result<This(A,BFL)> { typedef bf::cons<bool,typename boost::remove_reference<BFL>::type> type; };
 template<typename A, typename BFL> bf::cons<bool,BFL> operator ()(A &a, BFL const & bfl) const { return bf::cons<bool,BFL>(a.b(),bfl); }
};


int main() { 

 my_matrix_valued_function f1(1),f2(10);
 triqs::arrays::matrix<double> A(2,2); A()=0; A(0,0) = 2; A(1,1) = 20;

 TEST( (f1 - f2));
 TEST(f1(1));
 
 typedef struct triqs::arrays::matrix_expr<boost::proto::exprns_::basic_expr<boost::proto::tagns_::tag::terminal, boost::proto::argsns_::term<double>, 0l> > EE;
 typedef EE::proto_derived_expr yy;
 
 TEST( (2.0* f1 ) (0));
 TEST( (2* f1 + f2) (1));
 TEST( tqa::make_matrix ( ( 2*f1  +f2 ) (1)) );
 TEST( tqa::make_matrix  ( ( 2*f1  +f2 ) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  - f2 ) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  + 8*f2 ) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  - f2  - f2) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  + 8*f2 + f2 ) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  - f2  - 3*f2) (1)) );
 TEST( tqa::make_matrix ( ( 2*f1  + 8*f2 + f2/2 ) (1)) );
 
 // BUggy for the moment because of matmul to be generalized/...
 TEST( tqa::make_matrix ( ( A*f1  +f2 ) (1)) );

 print ( f1);
 print ( f1 - f2);

 //auto e = 2*f1  + 8*f2;
 //auto l = accumulator_transform<myaccu>()(e, bf::nil() );
 //std::cerr  << triqs::utility::typeid_name(l) <<std::endl;

};
