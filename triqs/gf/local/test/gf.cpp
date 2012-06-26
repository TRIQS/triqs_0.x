//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

#include <triqs/gf/local/gf.hpp>

//using namespace triqs::gf::local;
using namespace triqs::gf;
namespace tql= triqs::lazy;
//namespace tqa= triqs::arrays;
using tqa::range;

#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

std::complex<double> mul2 (std::complex<double> x) { return x*x;}

int main() {

 typedef local::gf<meshes::matsubara_freq> Gf_type;
 typedef local::gf_view<meshes::matsubara_freq> Gf_view_type;
 domains::freq_infty inf;

 Gf_type G1; // empty
 TEST( G1( 0) ) ;

 double beta =1;
 Gf_type G(2,2, meshes::matsubara_freq(beta,Fermion));
 Gf_type Gc(2,2, meshes::matsubara_freq(beta,Fermion));
 Gf_type G3(2,2, meshes::matsubara_freq(beta,Fermion));

 Gf_view_type Gv =G;
 TEST( G( 0) ) ;

 Gf_view_type Gv2 = slice(G,range(0,1),range(0,1));
 TEST( Gv2( 0) ) ;
 Gv2(0) = 10;
 TEST( Gv2( 0) ) ;
 TEST( G( 0) ) ;

 triqs::lazy::placeholder<0> om_;

 TEST( G(om_) ) ;
 TEST( tql::eval(G(om_), om_=0) ) ;

 TEST( Gv(om_) ) ;
 TEST( tql::eval(Gv(om_), om_=0) ) ;

 std::cout  <<"-------------lazy assign ------------------"<<std::endl;
 
 Gv(om_) = (0.2 + om_ + 2.1);
 TEST(G(0));
 TEST(G(inf));

 std::cout  <<"-------------lazy assign ------------------"<<std::endl;

 G(om_) = 1/(om_ + 2.3);

 TEST(G(0));
 TEST(G(inf));
 TEST(inverse(G(inf)));

 std::cout  <<"-------------------------------------"<<std::endl;

#define SUITE
#ifdef SUITE

 TEST( Gv(om_) ) ;
 TEST( tql::eval(Gv(om_), om_=0) ) ;

 // tail 
 BOOST_AUTO( t, G(inf));
 //local::gf<meshes::tail> t2 = t + 2.4;

 TEST(t.order_min()); 
 TEST( t( 0) ) ;

 TEST( Gv2(inf)( 0) ) ;

 // copy 
 Gc = G;
 TEST( G( 0) ) ;
 TEST( Gc( 0) ) ;


 // operations on gf
 G3 = G + Gc;

//#define ALL_TEST
#ifdef ALL_TEST
 for (int u=0; u<10; ++u) { 
  TEST( (G + 2.0* Gc)( u) ) ;
  TEST( (8.0*G + 2.0* Gc)( u) ) ;
  TEST( (8.0*G  - 2.0* Gc)( u) ) ;
  TEST( (G - Gc)( u) ) ;
  TEST( (G - 2.0* Gc)( u) ) ;
  TEST( (G * Gc)( u) ) ;
 }
#endif
 TEST( G( 0) ) ;
 TEST(G(inf)(0));

 TEST( ( G(inf) + G(inf) )  (0));
 TEST( ( G(inf) * G(inf) )  (0));
 TEST( (G + Gc)( inf) ) ;

 //TEST( (G + 2.3)(0));
 TEST( (t + 2.3)(0));

 TEST( t(-1));
 //TEST( (2 * inverse(t))(0));

 tqa::array<double,9> A(1,2,3,4,5,6,7,8,9);
 A()=0;
 //auto x = local::impl::gf_impl<triqs::gf::meshes::matsubara_freq, true>::wrap_infty (G.tail_view()) + 2.0;

 // test hdf5 
 H5::H5File file("ess_gf.h5", H5F_ACC_TRUNC );
 h5_write(file, "g", G);
#endif
}
