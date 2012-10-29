//#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK

#include <triqs/gf/lattice/glattice.hpp>

using namespace triqs::gf;
namespace tql= triqs::lazy;
using tqa::range;

#define TEST(X) std::cout << BOOST_PP_STRINGIZE((X)) << " ---> "<< (X) <<std::endl<<std::endl;

std::complex<double> mul2 (std::complex<double> x) { return x*x;}

int main() {

 typedef lattice::glattice<meshes::BZ_mesh,local::gf<meshes::matsubara_freq> > Gf_type;


 Gf_type G1; // empty
 TEST( G1( 0,0) ) ;

}
