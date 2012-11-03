//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

#include <triqs/gf/matsubara_freq.hpp> 
#include <triqs/gf/matsubara_time.hpp> 
#include <triqs/gf/block.hpp> 

namespace tql= triqs::clef;
namespace tqa= triqs::arrays;
using tqa::range;
using triqs::arrays::make_shape;
using triqs::gf::gf;
using triqs::gf::gf_view;
using triqs::gf::block;
using triqs::gf::Fermion;
using triqs::gf::matsubara_freq;
using triqs::gf::matsubara_time;

#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

int main() {

 double beta =1;
 auto G1 = matsubara_freq::make_gf (beta, Fermion, make_shape(2,2));
 auto G2 = matsubara_freq::make_gf (beta, Fermion, make_shape(2,2));
 auto G3 = matsubara_freq::make_gf (beta, Fermion, make_shape(2,2));

 std::vector<gf<matsubara_freq> >  V = { G1,G2,G3};
 std::vector<gf_view<matsubara_freq> >  Vv = { G1,G2,G3};

 std::cout <<" Building gf_view of view"<< std::endl ;
 auto GF_v = triqs::gf::block<matsubara_freq>::make_gf_view (Vv);

 std::cout <<" Building gf_view of gf"<< std::endl ;
 auto GF =  triqs::gf::block<matsubara_freq>::make_gf_view (V); //{G1,G2,G3});
 //auto GF = triqs::gf::block<matsubara_freq>::make_gf_view ( std::vector<gf_view<matsubara_freq> > {G1,G2,G3});

 auto  g0 = GF(0);
 auto  g0v = GF_v(0)();

 auto Gv = g0();

  Gv(0) = 20;
 TEST( Gv( 0) ) ;
 TEST( G1( 0) ) ;
 Gv(0) = 0;
  
 g0v(0) = 3.2;
  
// Vv[0](0) = -2.1;
 TEST( Gv( 0) ) ;
 TEST( G1( 0) ) ;


}
