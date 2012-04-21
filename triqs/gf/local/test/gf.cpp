#include <triqs/gf/local/gf.hpp>

//using namespace triqs::gf::local;
using namespace triqs::gf;
namespace tql= triqs::lazy;
namespace tqa= triqs::arrays;
#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

int main() {

 std::vector<std::vector<std::string> > ind(1); 
 ind[0].push_back("1");ind[0].push_back("2"); ind.push_back(ind[0]);

 typedef local::gf<meshes::matsubara_freq> Gf_type;
 typedef local::gf_view<meshes::matsubara_freq> Gf_view_type;

 Gf_type G1; // empty
 TEST( G1( 0) ) ;

 Gf_type G(2,2, meshes::matsubara_freq(Fermion), ind);
 Gf_type Gc(2,2, meshes::matsubara_freq(Fermion), ind);

 Gf_view_type Gv =G;
 TEST( G( 0) ) ;

 Gf_view_type Gv2 = G.slice(0,0);
 TEST( Gv2( 0) ) ;
 Gv2(0) = 10;
 TEST( Gv2( 0) ) ;
 TEST( G( 0) ) ;

 triqs::lazy::placeholder<0> om_;

 TEST( G(om_) ) ;
 TEST( tql::eval(G(om_), om_=0) ) ;

 TEST( Gv(om_) ) ;
 TEST( tql::eval(Gv(om_), om_=0) ) ;


  tqa::matrix<double> Id; Id() = 1;
  G(om_) = (om_ + 2.3);
  //G(om_) = 1/(om_ + 2.3);
 //G(om_) = Id* (1/(om_ + 2.3) );

 TEST( Gv(om_) ) ;
 TEST( tql::eval(Gv(om_), om_=0) ) ;

 // tail 
 BOOST_AUTO( t, G.tail_view());
 TEST( t( 0) ) ;

 TEST( Gv2.tail_view()( 0) ) ;

 // copy 
 Gc = G;
 TEST( G( 0) ) ;
 TEST( Gc( 0) ) ;

 TEST( (G + 2.0* Gc)( 0) ) ;
 TEST( (G * Gc)( 0) ) ;
 
 domains::infty inf;

 //meshes::matsubara_freq d = 
 TEST(G(inf)(0));

 TEST( ( G(inf) + G(inf) )  (0));
 TEST( ( G(inf) * G(inf) )  (0));
 TEST( (G + Gc)( inf) ) ;

 //auto x = local::impl::gf_impl<triqs::gf::meshes::matsubara_freq, true>::wrap_infty (G.tail_view()) + 2.0;
}
