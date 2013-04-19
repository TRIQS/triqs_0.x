//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

#include <triqs/gf/imfreq.hpp> 
#include <triqs/gf/imtime.hpp> 
#include <triqs/gf/block.hpp> 

namespace tql= triqs::clef;
namespace tqa= triqs::arrays;
using tqa::range;
using triqs::arrays::make_shape;
using triqs::gf::gf;
using triqs::gf::gf_view;
using triqs::gf::block;
using triqs::gf::Fermion;
using triqs::gf::imfreq;
using triqs::gf::imtime;
using triqs::gf::make_gf;
using triqs::gf::make_gf_view;

#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

int main() {

 double beta =1;
 
 auto G1 = make_gf<imfreq> (beta, Fermion, make_shape(2,2));
 auto G2 = make_gf<imfreq> (beta, Fermion, make_shape(2,2));
 auto G3 = make_gf<imfreq> (beta, Fermion, make_shape(2,2));

 std::vector<gf<imfreq> >  V ;
 V.push_back(G1); V.push_back(G2); V.push_back(G3); 
 std::vector<gf_view<imfreq> >  Vv; // = { G1,G2,G3};
 Vv.push_back(G1); Vv.push_back(G2); Vv.push_back(G3); 

 std::cout <<" Building gf_view of view"<< std::endl ;
 auto GF_v = make_gf_view<block<imfreq>> (Vv);

 std::cout <<" Building gf_view of gf"<< std::endl ;
 auto GF =  make_gf_view<block<imfreq>> (V); //{G1,G2,G3});
 //auto GF = make_gf_view<block<imfreq>> ( std::vector<gf_view<imfreq> > {G1,G2,G3});

 std::cout  << "Number of blocks " << GF.mesh().size()<<std::endl ;
 auto g0 = GF[0];
 auto g0v = GF_v[0]();

 auto Gv = g0();

 Gv.on_mesh(0) = 20;
 TEST( Gv( 0) ) ;
 TEST( G1( 0) ) ;
 Gv.on_mesh(0) = 0;

 g0v.on_mesh(0) = 3.2;
 //g0v[0]= 3.2;

 // Vv[0](0) = -2.1;
 TEST( Gv( 0) ) ;
 TEST( G1( 0) ) ;
 //TEST( GF_v(0)( 0) ) ;
 //TEST( GF_v[0]( 0) ) ;

 // bug fixed for this
 gf< block< imfreq> > G9;
 G9 = make_gf<block<imfreq>> (2, make_gf<imfreq>(beta, Fermion, make_shape(2,2)));

 // Operation
 g0.on_mesh(0) = 3.2;
 TEST( GF[0](0) ) ;
 GF = GF/2; 
 TEST( GF[0](0) ) ;
 //TEST( g0("3.2") ) ;
 //TEST( GF(0)(0) ) ;

}
