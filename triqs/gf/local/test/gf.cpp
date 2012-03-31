#include <triqs/gf/local/gf.hpp>

//using namespace triqs::gf::local;
using namespace triqs::gf;
namespace tql= triqs::lazy;
#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

int main() {

 std::vector<std::string> ind; ind.push_back("1");ind.push_back("2");
 typedef local::gf<meshes::matsubara_freq,false> Gf_type;
 typedef local::gf<meshes::matsubara_freq,true> Gf_view_type;
 
 Gf_type G1; // empty
 TEST( G1( 0) ) ;

 Gf_type G(2,2, meshes::matsubara_freq(Fermion), ind,ind);
 
 // should not compile and does not
 //Gf_view_type Ge(2,2, meshes::matsubara_freq(Fermion), ind,ind);

 Gf_view_type Gv =G;
 TEST( G( 0) ) ;

 Gf_view_type Gv2 = G.slice(0,0);
 TEST( Gv2( 0) ) ;

 triqs::lazy::placeholder<0> om_;

 TEST( G(om_) ) ;
 TEST( tql::eval(G(om_), om_=0) ) ;

 TEST( Gv(om_) ) ;
 TEST( tql::eval(Gv(om_), om_=0) ) ;

 G(om_) = 1/(om_ + 2.3);

  TEST( Gv(om_) ) ;
 TEST( tql::eval(Gv(om_), om_=0) ) ;

 // tail 
 BOOST_AUTO( t, G.tail_view());
 TEST( t( 0) ) ;

 TEST( Gv2.tail_view()( 0) ) ;
  
}
